#Sys.setenv(KBPIPEOUTDIR='/home/kevbrick/helix/RDCO/kevbrick/GL_Sorting_Paper/dataReRun2/output/')
#Sys.setenv(KBPIPEWORKDIR='/home/kevbrick/RDCO/kevbrick/GL_Sorting_Paper/')

## Load custom functions & H3 data
workdir <- Sys.getenv('KBPIPEWORKDIR')
outdir  <- Sys.getenv('KBPIPEOUTDIR')
datadir <- paste0(outdir,'/Rtables/')

if (exists('testingKBGLPIPE')){
  imgdir = '/home/kevbrick/testImg/'
}else{
  #imgdir  <- paste0(outdir,'/figs/')
  imgdir  <- './'
}

## Load custom functions & H3 data
#source(paste0(workdir,'accessoryFiles/scripts/R/sortingPaper_functions.R'))
source(paste0(workdir,'accessoryFiles/scripts/R/genericFunctionsLamEtAl.R'))
source(paste0(workdir,'accessoryFiles/scripts/R/sortingPaper_plotExpressionForSingleGene_v_Bortvin.R'))

library(ggplot2)
library(reshape2)
library(plyr)
library(factoextra)

eData   <- read.table(paste0(datadir,'tssTable.forR.tab'),header=TRUE)

eData$type                                           <- 'Other'
eData$type[eData$TSS & !eData$HS]                    <- 'TSS'
eData$type[eData$TSS & eData$HS]                     <- 'HSTSS'
eData$type[eData$HS & !eData$HSPrdm9KO & !eData$TSS] <- 'HS'

#Merge EP and LP H3K4me3 data
eData$PA.R1.H3K4me3.NCISTSS <- (eData$EP.R1.H3K4me3.NCISTSS + 
                                eData$LP.R1.H3K4me3.NCISTSS)/2 

#eData$PA.R1.H3K4me3.KM      <- (eData$EP.R1.H3K4me3.KM +
#                                eData$LP.R1.H3K4me3.KM)/2

##Assign column names & labels
expLabel <- c('L','Z','P','D')
expNames <- paste0('e',c('LE','ZY','PA','DI'))

stageLabel <- c("L","Z","P","D")
#stageNames <- c("LE.R1.H3K4me3.KM",
#                "ZY.R1.H3K4me3.KM",
#                "PA.R1.H3K4me3.KM",
#                "DI.R1.H3K4me3.KM")

stageNames <- c("LE.R1.H3K4me3.NCISTSS",
                "ZY.R1.H3K4me3.NCISTSS",
                "PA.R1.H3K4me3.NCISTSS",
                "DI.R1.H3K4me3.NCISTSS")


## Subset down to only columns of interest
dTSS <- eData[eData$type == 'TSS',
              c("cs","from","to",
                "type","name","numExp",
                stageNames,expNames)]

## Remove miRNA, rRNA (no longer an issue, but retain)
dTSS <- dTSS[!(grepl('Mir',dTSS$name)),]
dTSS <- dTSS[!(grepl('n-R5s',dTSS$name)),]
dTSS <- dTSS[!(grepl('mt-Rnr',dTSS$name)),]

## Set ChIP-Seq <0 to 0
for (n in stageNames){
  dTSS[dTSS[,n]<=0,n] <- 0
}

## Require > 0.2 TPM expression
dTSS$TPMOK    <- apply(dTSS[,expNames], 1, max) >= 0.2

## No zero/neg vals for expression
dTSS$expZero  <- apply(dTSS[,expNames],1,function(x){(sum(x<=0)<1)})

## No zero/neg vals for H3K4me3
dTSS$nZero    <- apply(dTSS[,stageNames],1,function(x){(sum(x<=0)<1)})

## 1.2x fold change
dTSS$FC  <- apply(dTSS[,expNames],1,function(x){(max(x)/min(x[x>0]))>1.25})

## Apply filters
dTSS <- dTSS[dTSS$FC,]
dTSS <- dTSS[dTSS$TPMOK,]
dTSS <- dTSS[dTSS$expZero,]
dTSS <- dTSS[dTSS$nZero,]

## For both H3K4me3 AND RNA-Seq
## Normalize by mean - s.d. then set range from 0->1

dTSS[,expNames]   <- standardizeMNSD(dTSS[,expNames])
dTSS[,stageNames] <- standardizeMNSD(dTSS[,stageNames])

useOnlySingleTranscriptGenes <- TRUE
if (useOnlySingleTranscriptGenes){
  dTSS <- dTSS[dTSS$numExp <= 2,]
}

## Remove NANs and Infs
for (n in stageNames){
  dTSS[is.nan(dTSS[,n]),n] <- 0
  dTSS[is.infinite(dTSS[,n]),n] <- 0
}

dTSS$noNAN <- !is.nan(apply(dTSS[,expNames],1,max))
dTSS       <- dTSS[dTSS$noNAN,]

## Get "optimal" number of clusters
getBestClusterSize <- TRUE;

if (getBestClusterSize){
  fvH3 <- fviz_nbclust(dTSS[,stageNames],
                       kmeans,
                       method = "gap_stat",
                       nboot = 100,
                       k.max = 18)

  K_H3  <- fvH3$data$gap - c(0,fvH3$data$gap[1:(length(fvH3$data$gap)-1)])
  nK_H3 <- min(which (K_H3<0))-1

}else{
  print (" *** WARNING *** Best cluster detection DISABLED !!! This is only for rapid testing ... turn it BACK ON !!!\n")
  bestK <- 5
}

if (exists('fvH3')){
  optimalK <- which(!(fvH3$data$gap > c(0,fvH3$data$gap[1:(length(fvH3$data$gap)-1)])))[1]-1
}else{
  optimalK <- 5
}

for (nK in optimalK){
  ## Do k-means clustering
  k5 <- kmeans(dTSS[,stageNames],centers = nK,iter.max = 50)
  
  ## Add to DF
  dTSS$k5       <- k5$cluster
  dTSS$cluster  <- k5$cluster
  
  #dTSS <- reOrderClusters(dTSS,expNames)
  dTSS <- reOrderClusters(dTSS,stageNames)
  
  ## Calculate cluster medians
  a <- aggregate(x = dTSS,
                 by = list(dTSS$k5),
                 FUN=function(x){median(as.numeric(x),na.rm=TRUE)})
  
  aSD <- aggregate(x = dTSS,
                 by = list(dTSS$k5),
                 FUN=function(x){sd(as.numeric(x),na.rm=TRUE)/sqrt(length(x))})
  
  ## Set cluster names with N values
  names(a)[which(names(a)=='k5')]   <- 'statK'
  names(aSD)[which(names(a)=='k5')] <- 'statK'
  names(a)[1]                       <- 'k5'
  names(aSD)[1]                     <- 'k5'
  
  a$cluster                         <- 'X'
  aSD$cluster                       <- 'X'
  clusterOrder                      <- matrix('Unknown',nK)
  
  for (t in 1:nK){
    clusterName                       <- paste0("(",t,") N=",format(sum(dTSS$k5 == t), big.mark=","))
    clusterOrder[t]                   <- clusterName
    dTSS$cluster[which(dTSS$k5==t)]   <- clusterName
    a$cluster[which(a$k5==t)]         <- clusterName
    aSD$cluster[which(aSD$k5==t)]     <- clusterName
  }
  
  ## Set order of clusters (by peak time)
  dTSS$cluster  <- factor(dTSS$cluster,levels=clusterOrder)
  a$cluster     <- factor(a$cluster,levels=clusterOrder)
  aSD$cluster   <- factor(aSD$cluster,levels=clusterOrder)
  
  mH3K4 <- reshape2:::melt.data.frame(a,id.vars = c("k5","cluster"),
                                      measure.vars = stageNames )
  mExpr <- reshape2:::melt.data.frame(a,id.vars = c("k5","cluster"),
                                      measure.vars = expNames )

  sdH3K4 <- reshape2:::melt.data.frame(aSD,id.vars = c("k5","cluster"),
                                       measure.vars = stageNames )
  sdExpr <- reshape2:::melt.data.frame(aSD,id.vars = c("k5","cluster"),
                                       measure.vars = expNames )
  
  mH3K4$sd <- sdH3K4$value
  mExpr$sd <- sdExpr$value
  
  ## number all genes
  dTSS$gene <- 1:dim(dTSS)[1]
  
  ## Do random permutations for correlation coefficient
  ccFigs <- makeCCFigs(dTSS,nK,100)
  
  ## add p-values to DF
  dTSSwithP <- plyr:::join(dTSS,ccFigs$pvals,by="newClusterNumber")
  
  allH3K4   <- reshape2:::melt.data.frame(dTSSwithP,id.vars = c("cluster","p","cc","gene"),measure.vars = stageNames)
  allExpr   <- reshape2:::melt.data.frame(dTSSwithP,id.vars = c("cluster","p","cc","gene"),measure.vars = expNames)

  ## Not used now: numRows <- ceiling(nK/6)
  numRows <- 1
  
  g1 <- ggplot(allH3K4[allH3K4$gene<2500,],aes(x=variable,y=value)) + 
    facet_wrap(~cluster,nrow=numRows) + 
    geom_line(aes(group=gene),
              color='grey50',alpha=.015) + 
    geom_line(data=mH3K4,size=0.2,
              aes(x=variable,y=value,
                  group=cluster),
              color='black') + 
    geom_point(data=mH3K4,size=0.9,
               aes(x=variable,y=value,
                   group=cluster),
               color='black',
               shape=21,
               fill='orange') + 
    theme(legend.position='none',
          strip.background = element_blank(),
          strip.text = element_blank()) + 
    xlab('Stage') + ylab('H3K4me3\nsignal at\nTSS (AU)') + 
    scale_x_discrete(breaks=stageNames,labels=stageLabel) 

  g2 <- ggplot(allExpr[allExpr$gene<2500,],aes(x=variable,y=value)) + 
    facet_wrap(~cluster,nrow=numRows) +
    geom_line(aes(group=gene),
              color='grey50',alpha=.01) + 
    geom_line(data=mExpr,size=0.2,
              aes(x=variable,y=value,
                  group=cluster),
              color='black') +
    geom_point(data=mExpr,size=0.9,
               aes(x=variable,y=value,
                   group=cluster),
               color='black',
               shape=21,
               fill='orange') +  
    theme(legend.position='none',
          strip.background = element_blank(),
          strip.text = element_blank()) + 
    xlab('Days post-partum (dpp)') + 
    ylab('Transcript\nExpression\n(RNA-Seq; AU)') + 
    scale_x_discrete(breaks=expNames,labels=expLabel) + 
    scale_y_continuous(breaks=c(-1,0,1)) + 
    geom_text(x=-Inf,y=-Inf,
              aes(label=paste(paste0(" ",cluster),paste0('R = ',round(cc,2)),sep='; ')),
              vjust=-1, hjust=0, size=7*5/14,
              check_overlap = TRUE)
}

noMargin <- theme(plot.margin = unit(c(0,0,0,0),'cm'))
noX      <- theme(axis.text.x = element_blank())

gK4vExp <- ggarrange(g1 + noMargin + noX + xlab (''),
                g2 + noMargin + xlab('MPI stage'),
                nrow=2, ncol=1,
                align='v')

### Get CC for ALL genes and p-val
nIter <- 10000
lCC <- getH3vExpCC(dTSS,nIter)

gSim <- ggplot(as.data.frame(lCC$geneShuf),aes(x=V1)) + 
  geom_density(color=NA,alpha=.4,fill='steelblue') + 
  geom_vline(xintercept=lCC$cc,color='magenta') + 
  xlab('Correlation coefficient') + 
  ylab('Density') + 
  annotate(geom='text',label=paste0("All transcripts R2 = ",round(lCC$cc,2),"   "),
           x=lCC$cc,
           y=Inf,
           vjust=1,
           color='magenta',
           hjust=1,
           size=8*5/14) + 
  annotate(geom='text',
           label=paste0(comma(nIter)," bootstraps\n","mean(R2) = ",round(mean(lCC$geneShuf),2),"\n\n"),
           x=max(lCC$geneShuf)*2,
           y=0,
           vjust=0,
           color='steelblue',
           hjust=0.5,
           size=8*5/14)

png(getIMGname(fname = 'H3K4me3_v_Expression_AllGenes_v_Bootstraps',
               type = 'PNG',
               saveLocation = imgdir),
    width=3,
    height=2,
    units='in',
    res=400)
print(gSim)
dev.off()

png(getIMGname(fname = 'H3K4me3_v_Expression_AllGenes_v_Bootstraps',
               type = 'PDF',
               saveLocation = imgdir),
    width=3,
    height=2)
print(gSim)
dev.off()


###############################

