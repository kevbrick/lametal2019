Sys.setenv(KBPIPEOUTDIR='/data/RDCO/kevbrick/GL_Sorting_Paper/finalPipe/output/')
Sys.setenv(KBPIPEWORKDIR='/data/RDCO/kevbrick/GL_Sorting_Paper/')

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

print(workdir)

## Load custom functions & H3 data
source(paste0(workdir,'accessoryFiles/scripts/R/genericFunctionsLamEtAl.R'))

library(ggplot2)
library(ggrepel)
library(scales)
library(factoextra)

theme7point()

## Begin
eData <- read.table(paste0(datadir,'/abk4meTable.forR.tab'),header=TRUE)

countOLs <- function (x){
  
  aggX2            <- aggregate(x$type,by=list(x$cluster,x$type),FUN=function(x){length(x)})
  names(aggX2)     <- c('cluster','newType','N')
  
  aggTotC          <- aggregate(x$type,by=list(x$cluster),FUN=function(x){length(x)})
  names(aggTotC)   <- c('cluster','cTot')
  
  aggTotT          <- aggregate(x$type,by=list(x$type),FUN=function(x){length(x)})
  names(aggTotT)   <- c('newType','tTot')
  
  aggBoth          <- join(aggX2,aggTotC,"cluster")
  aggX3            <- join(aggBoth,aggTotT,"newType")
  
  aggX3$clustPC    <- aggX3$N/aggX3$cTot
  aggX3$clustPCLbl <- paste0(' ',round(aggX3$clustPC*100,0),'%')
  
  aggX3$typePC     <- aggX3$N/aggX3$tTot
  aggX3$typePCLbl  <- paste0(' ',round(aggX3$typePC*100,0),'%')
  
  return(aggX3)
}

makeTables <- function(x){
  #############################################
  #dfDefault         <- x[x$type !=%in% c('KO','TSSKO'),]
  dfDefault         <- x[x$type != "Hotspot" & x$type != "HS & TSS",]
  dfDefault$type    <- factor(dfDefault$type,levels=c("TSS","Other","Non-HS","Default"))
  dfDefault$type    <- 'Non-HS'
  dfDefault$type[dfDefault$HSPrdm9KO] <- 'Default'
  
  dfAll             <- x[x$type %in% c('Other','TSS'),]
  dfAll$type     <- 'Non-HS'
  
  aggDF      <- countOLs(x)
  aggDefault <- countOLs(dfDefault)
  aggAll     <- countOLs(dfAll)
  
  aggJustDef <- rbind(aggDF,aggDefault,aggAll)
  
  noleg  <- theme(legend.position='none')
  noMarg <- theme(plot.margin=unit(c(0,0,0,0),'cm'))
  myColorScale <- scale_fill_gradient2(low='white',mid='tomato',high='firebrick',
                                       midpoint=.5,
                                       limits=c(0,1)) 
  
  ggClust <- ggplot(aggDF,
                    aes(x=cluster, y=newType,
                        fill=clustPC)) +
    geom_tile() + 
    ylab('') + xlab('') +
    geom_text(aes(label=clustPCLbl,size=clustPC)) + 
    scale_size(range = c(5*5/14, 7*5/14), guide = F) +
    myColorScale +  
    theme(axis.text.x=element_blank()) +
    coord_cartesian(ylim=c(0.5,4.5),expand=FALSE) + 
    ggtitle('Cluster composition') + noleg + noMarg
  
  ggType <- ggplot(aggDF,
                   aes(x=cluster, y=newType,
                       fill=typePC)) +
    geom_tile() + 
    geom_text(aes(label=typePCLbl,
                  size=typePC)) + 
    coord_cartesian(ylim=c(0.5,4.5),expand=FALSE) + 
    scale_size(range = c(5*5/14, 7*5/14), guide = F) +
    myColorScale +  
    ggtitle('Cluster distribution of H3K4me3 peak sub-types') +
    noleg + noMarg + 
    xlab('') + ylab('') + theme(axis.text.x=element_blank())
  
  ggTypeDEF <- ggplot(aggDefault,
                      aes(x=cluster, y=newType,
                          fill=typePC)) +
    geom_tile() + 
    geom_text(aes(label=typePCLbl,size=typePC)) + 
    scale_size(range = c(5*5/14, 7*5/14), guide = F) +
    coord_cartesian(ylim=c(0.5,2.5),expand=FALSE) + 
    myColorScale +  
    noleg + noMarg + ylab('') + 
    ggtitle('Cluster distribution of sites used as "default" hotspots') +
    scale_y_discrete(breaks=c('Non-HS','Default'),labels=c('Not used','Default HS'))
  
  gReturn <- ggarrange(ggClust,
                       ggType ,
                       ggTypeDEF,
                       heights    = c(1.6,1.6,1.3),
                       labels     = c('c','d','e'),
                       vjust      = 1,
                       hjust      = 0,
                       align      = 'v',
                       font.label = list(size = 8,face='bold'),
                       ncol       = 1,
                       nrow       = 3)
  
  return(list(gAll=gReturn,
              gTop =ggClust,
              gMid =ggType,
              gEnd = ggTypeDEF))
  
}

stageLabel <- c("LE","ZY","EP","LP","DI")

stageNames <- paste0(c("LE","ZY","EP","LP","DI"),".R2.abK4me.NCISTSS")

eData$type <- 'Other'
eData$type[eData$TSS & !eData$HS] <- 'TSS'
eData$type[!eData$TSS & eData$HS] <- 'Hotspot'
eData$type[eData$TSS & eData$HS] <- 'HS & TSS'

eData$type <- factor(eData$type, 
                     levels=c('Hotspot','HS & TSS','TSS','Other'))

eData$csType <- 'other'
eData$csType[eData$cs == "chrX"] <- 'chrX'
eData$csType[eData$cs %in% paste0("chr",1:19)] <- 'autosome'

dTypeCount <- as.data.frame(table(eData$type))
names(dTypeCount) <- c('type','count')
nTotal <- sum(dTypeCount$count)

dTypeCount$pc    <- dTypeCount$count/nTotal*100
dTypeCount$yPos  <- 100 - (dTypeCount$pc/2 + c(0, cumsum(dTypeCount$pc)[-length(dTypeCount$pc)]))
dTypeCount$myLbl <- paste0(dTypeCount$type,"; ",percent(dTypeCount$pc/100),"\nN = ",comma(dTypeCount$count))

## Panel A: pie chart
gPie <- ggplot(dTypeCount,
               aes(x="",y=pc,fill=type)) + 
  geom_bar(width=1, stat="identity", color='black',lwd=.2) + 
  geom_label_repel(aes(y = yPos, 
                       label = myLbl), 
                   size=7*5/14) + 
  coord_polar('y') + 
  theme(axis.text=element_blank(),
        axis.ticks = element_blank(), 
        axis.line = element_blank(),
        legend.position='none') + 
  scale_fill_manual(values=c('lightgreen','darkseagreen3','grey60','khaki1')) +
  xlab('') + ylab('')

## Allow "Default" as a type 
eData$type <- factor(eData$type, 
                     levels=c('Hotspot','HS & TSS','TSS','Default','Other'))

dfClust <- eData[eData$type %in% c('TSS','Hotspot','HS & TSS','Other'),]

## Number the genes
dfClust$id <- 1:length(dfClust$cs)

dfClust[,stageNames] <- standardizeMNSD(dfClust[,stageNames])

dfClust <- dfClust[notNanRows(dfClust[,stageNames]),]

getBestClusterSize <- TRUE;
if (getBestClusterSize == TRUE){
  fvAll <- fviz_nbclust(dfClust[,stageNames], 
                        kmeans, 
                        method = "gap_stat", 
                        k.max = 12, nboot = 50)
  
  K_H3  <- fvAll$data$gap - c(0,fvAll$data$gap[1:(length(fvAll$data$gap)-1)]) 
  nK_H3 <- min(which (K_H3<0))-1
}

if (exists('fvAll')){
  nK <- which(!(fvAll$data$gap > c(0,fvAll$data$gap[1:(length(fvAll$data$gap)-1)])))[1]-1
}else{
  nK <- 5
}

## Do k-means clustering
k5 <- kmeans(dfClust[,stageNames],centers = nK,iter.max = 50)

dfClust$k5       <- k5$cluster

## Renumber clusters to retain pseudo-temporal order
dfClust <- reOrderClusters(dfClust,stageNames)

remove(clusterOrder)
remove(clusterName)

clusterOrder <- rep('',nK)

## Rename clusters with N
for (t in 1:nK){
  clusterName                           <- paste0("(",t,")"," N=",format(sum(dfClust$k5 == t), big.mark=","))
  clusterOrder[t]                       <- clusterName
  dfClust$cluster[which(dfClust$k5==t)] <- clusterName
}

## Melt DF to allow for plotting of all lines
mdfClust <- reshape2:::melt.data.frame(dfClust,
                                       id.vars=c('id','cluster'),
                                       measure.vars=stageNames,
                                       variable.name = 'stage',
                                       value.name = "strength")

mdfClust$stage <- factor(mdfClust$stage, levels=c(stageLabel,stageNames))

for (nS in 1:length(stageNames)){
  mdfClust$stage[mdfClust$stage==stageNames[nS]] <- stageLabel[nS]
}

## Calculate cluster medians
meanClust <- aggregate(mdfClust$strength,
                       by = list(mdfClust$cluster,
                                 mdfClust$stage),
                       FUN=function(x){
                         mean(as.numeric(x),na.rm=TRUE)
                       }
)

names(meanClust) <- c('cluster','stage','strength')

## Set order of clusters
dfClust$cluster   <- factor(dfClust$cluster,levels=clusterOrder)
mdfClust$cluster  <- factor(mdfClust$cluster,levels=clusterOrder)
meanClust$cluster <- factor(meanClust$cluster,levels=clusterOrder)

## Make lables DF
labels <- as.data.frame(table(meanClust$cluster))[1]
names(labels) <- 'cluster'
labels$xVal   <- 'LE'

## Panel B: Draw profiles for all genes
gProfiles <- ggplot(mdfClust,aes(x=stage,y=strength)) +
  facet_wrap(~cluster,ncol=nK) + 
  geom_hline(yintercept=0,color='grey30',lwd=0.2) +
  geom_line(data=mdfClust[mdfClust$id %in% sample(1:(dim(dfClust)[1]),10000),],
            aes(group=id), lwd=.2, alpha=0.03) + 
  geom_line(data=meanClust,
            aes(group=cluster),
            color='mediumpurple1',
            lwd=.3) + 
  geom_point(data=meanClust,
             size=1,
             shape=21,
             color='mediumpurple1',
             fill='mediumpurple3') + 
  theme(legend.position='none',
        strip.background=element_blank(),
        strip.text=element_blank()) + 
  geom_text(data=labels, 
            y=-Inf, 
            hjust=0, vjust=-0.5, 
            size=7*5/14, 
            aes(x=xVal, 
                label=paste0(" ",cluster)),
            check_overlap = TRUE,
            color='black') + 
  xlab('Stage') + ylab('H3K4me2/3')  

## Do random permutations for correlation coefficient
#ccFigs <- makeCCFigs(dTSS,nK,5)

## make heatmaps
gOut <- makeTables(dfClust)
gOutA <- makeTables(dfClust[dfClust$csType == 'autosome',])

################################################3
panelA <- gPie + ggtitle('H3K4me2/3 peaks by type') + 
  theme(plot.title = element_text(hjust=0.5),
        plot.margin=unit(c(0,0,0,0),'cm'))

panelB <- gProfiles

panelC <- gOut$gAll

ggAC <- ggarrange(panelA,panelC,
                  ncol=2,nrow=1,
                  labels=c('a','c'),
                  vjust = 1,hjust=0,
                  widths=c(1,1.5),
                  font.label = list(size = 8,face='bold'))

ggABC <- ggarrange(ggAC,panelB,
                   ncol=1,nrow=2,
                   labels=c('','b'),
                   vjust = 1,hjust=0,
                   heights=c(1,.6),
                   font.label = list(size = 8,face='bold'))

png(getIMGname(fname = 'Figure_S7',
               type = 'PNG',
               saveLocation = imgdir),
    width=7,
    height=4,
    units='in',
    res=400)
print(ggABC)
dev.off()

pdf(getIMGname(fname = 'Figure_S7',
               type = 'PDF',
               saveLocation = imgdir),
    width=7,
    height=4)
print(ggABC)
dev.off()
