#Sys.setenv(KBPIPEOUTDIR='/home/kevbrick/helix/RDCO/kevbrick/GL_Sorting_Paper/dataReRun2/output/')
#Sys.setenv(KBPIPEWORKDIR='/home/kevbrick/helix/RDCO/kevbrick/GL_Sorting_Paper/')

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
source(paste0(workdir,'accessoryFiles/scripts/R/genericFunctionsLamEtAl.R'))
source(paste0(workdir,'accessoryFiles/scripts/R/sortingPaper_plotExpressionForSingleGene_v_Bortvin.R'))

library(ggplot2)
library(reshape2)
library(plyr)
library(purrr)

theme7point()

clusterCluster <- function(data2Clust,
                           eNames,
                           sNames,
                           xLabels = c('LE','ZY','PA','DI'),
                           originalCluster = NULL,
                           nGenes2Show = 6){
  
  if (is.null(originalCluster)){originalCluster <- data2Clust$nK[1]}
  
  for (i in 1:4){
    data2Clust[[paste0('k4',as.character(xLabels[i]))]]  <- data2Clust[[as.character(sNames[i])]] 
  }
  
  myTSS <- data2Clust
  data2Clust <- data2Clust[ , !(names(data2Clust) %in% c('nK','cluster'))]
  
  fvExp <- fviz_nbclust(data2Clust[,eNames], 
                        kmeans, 
                        method = "gap_stat", 
                        k.max = 16)
  
  K_exp <- fvExp$data$gap - c(0,fvExp$data$gap[1:(length(fvExp$data$gap)-1)]) 
  nK <- min(which (K_exp<0))-1

  k5 <- kmeans(data2Clust[,expNames],centers = nK,iter.max = 50)
  data2Clust$k5 <- k5$cluster

  data2Clust <- reOrderClusters(data2Clust,eNames)
  
  ## Calculate cluster medians
  a <- aggregate(x = data2Clust,
                 by = list(data2Clust$k5),
                 FUN=function(x){median(as.numeric(x),na.rm=TRUE)})
  
  aSD <- aggregate(x = data2Clust,
                   by = list(data2Clust$k5),
                   FUN=function(x){sd(as.numeric(x),na.rm=TRUE)/sqrt(length(x))})
  
  ## Set cluster names with N values
  names(a)[which(names(a)=='k5')] <- 'statK'
  names(aSD)[which(names(a)=='k5')] <- 'statK'
  names(a)[1] <- 'k5'
  names(aSD)[1] <- 'k5'
  
  a$cluster <- 'X'
  aSD$cluster <- 'X'
  clusterOrder <- matrix('Unknown',nK)
  
  for (t in 1:nK){
    clusterName                       <- paste0("N=",format(sum(data2Clust$k5 == t), big.mark=","))
    clusterOrder[t]                   <- clusterName
    data2Clust$cluster[which(data2Clust$k5==t)]   <- clusterName
    a$cluster[which(a$k5==t)]         <- clusterName
    aSD$cluster[which(aSD$k5==t)]     <- clusterName
  }
  
  data2Clust$cluster <- factor(data2Clust$cluster,levels=clusterOrder)
  
  ## Set order of clusters
  #data2Clust$cluster  <- factor(data2Clust$cluster,levels=clusterOrder)
  #a$cluster           <- factor(a$cluster,         levels=clusterOrder)
  #aSD$cluster         <- factor(aSD$cluster,       levels=clusterOrder)
  
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
  
  mExpr$type <- 'RNA-Seq'
  
  data2Clust$gene <- 1:dim(data2Clust)[1]
  
  ## Do random permutations for correlation coefficient
  ccFigs <- makeCCFigs(data2Clust,nK,1000)
  
  dwithP <- plyr:::join(data2Clust,ccFigs$pvals,by="newClusterNumber")
  allH3K4   <- reshape2:::melt.data.frame(dwithP,id.vars = c("cluster","p","cc","gene"),measure.vars = stageNames)
  allExpr   <- reshape2:::melt.data.frame(dwithP,id.vars = c("cluster","p","cc","gene"),measure.vars = expNames)
  
  allH3K4$variable <- factor(allH3K4$variable,levels=c(stageNames,xLabels))

  allH3K4$variable[allH3K4$variable == sNames[1]] <- 'LE'
  allH3K4$variable[allH3K4$variable == sNames[2]] <- 'ZY'
  allH3K4$variable[allH3K4$variable == sNames[3]] <- 'PA'
  allH3K4$variable[allH3K4$variable == sNames[4]] <- 'DI'

  allH3K4$type <- 'H3K4me3'
  
  allExpr$variable <- factor(allExpr$variable,levels=c(expNames,xLabels))
  
  allExpr$variable[allExpr$variable == 'eLE'] <- 'LE'
  allExpr$variable[allExpr$variable == 'eZY'] <- 'ZY'
  allExpr$variable[allExpr$variable == 'ePA'] <- 'PA'
  allExpr$variable[allExpr$variable == 'eDI'] <- 'DI'
  
  allExpr$type <- 'RNA-Seq'
  
  allBoth <- rbind(allExpr,allH3K4)
  ## Not used now
  #numRows <- ceiling(nK/6)
  numRows <- 1
  
  allBoth$title <- paste(paste0(" ",allBoth$cluster),
                         paste0('R = ',round(allBoth$cc,2)),
                         allBoth$p,
                         sep='; ')
  
  lblOrder <- matrix(nrow=1,ncol=nK)
  for (k in 1:nK){
    lblOrder[k] <- allBoth$title[allBoth$cluster == clusterOrder[k]][1]
  }
  
  allBoth$title <- factor(allBoth$title,levels=lblOrder)
  
  gInit <- ggplot(allBoth[allBoth$gene<2500,],
         aes(x=variable,y=value)) + 
    geom_line(aes(group=paste0(gene,type),color=type),
              alpha=.01) +
    theme(legend.position='none',
          legend.title=element_blank(),
          strip.background = element_blank()) + 
    xlab('') + 
    ylab('') + 
    scale_x_discrete(breaks=xLabels,labels=xLabels) + 
    scale_y_continuous(breaks=c(-1,0,1)) + facet_wrap(~type)
  
  gC2 <- ggplot(allBoth[allBoth$gene<2500,],
                 aes(x=variable,y=value)) + 
    facet_grid(type~title) +
    geom_line(aes(group=paste0(gene,type),color=type),
              alpha=.01) + 
    theme(legend.position='none',
          strip.background = element_blank(),
          strip.text.y = element_blank(),
          strip.text.x= element_text(size=7)) + 
    xlab('') + 
    ylab('') + 
    scale_x_discrete(breaks=xLabels,labels=xLabels) + 
    scale_y_continuous(breaks=c(-1,0,1)) 
  
  ## For Debugging
  if (any(names(data2Clust) == 'oldname')){
    data2Clust$name <- data2Clust$oldname  
  }
  
  data2Clust$oldname <- data2Clust$name
  
  ## Parse out gene name
  data2Clust$name <- map_chr(data2Clust$oldname, ~as.vector(strsplit(as.character(.x[1]), split = '[|]'))[[1]][6])
  
  meiGenes <- c('Spo11','Mei4','Smc3','Smc1',
                'Sycp3','Sycp1','Histh1t','Prdm9',
                'Tex11','Tex19','Rec114','Hormad1',
                'Dmc1','Hormad2','Brca1','Rad51',
                'Stag3','Rec8','Ankrd31','Ctcfl,Ctcflos',
                'Sycp2','Sycp2l','Iho1','Stra8','Mlh1','Mlh3',
                'Msh4','Msh2','Msh6','Cntd1','Gm960','Ccdc36',
                'Mei1','Meiob','Spata22','Wdr61','Mre11','Nbs1','Rad50',
                'Rad54')
  
  mG       <- data.frame(name=meiGenes)
  
  remove(gAll)
  for (nCluster in 1:nK){
    remove(pAll)
    
    dC <- data2Clust[data2Clust$k5 == nCluster,]
    
    dCM <- plyr:::join(dC,mG,"name",'inner')
    
    nOther <- nGenes2Show-dim(dCM)[1]
    
    if (nOther > 0){
      nameGene <- c(as.character(dCM$name),as.character(sample(data2Clust$name[data2Clust$k5 == nCluster],nOther)))
    }else{
      nameGene <- as.character(sample(dCM$name,nGenes2Show,replace=FALSE))
    }
    
    noLeg <- theme(legend.position='none')
    noMarg <- theme(plot.margin=unit(c(0,0,0,0),'cm'))
    
    pAll <- plotGeneExpNFMulti(nameGene,data2Clust,noXlab = TRUE,noYlab = TRUE,tall = TRUE) + 
      noLeg + noMarg + 
      scale_y_continuous(breaks=c(-1,0,1)) 
    
    if (exists('gAll')){
      h1   <- nCluster-1
      h2   <- nCluster-h1
      gX  <- gAll
      gAll <- ggarrange(gX,pAll,ncol=2,nrow=1,widths=c(h1,h2))
    }else{
      gAll <- pAll
    }
  }  
  
  yTitle <- textGrob("Relative signal", gp=gpar(fontsize=7), rot=90)
  xTitle <- textGrob("MPI stage", gp=gpar(fontsize=7))
  
  marginOK <- theme(plot.margin = unit(c(0.2,0,0,0),units = 'cm'))
  
  gAInit   <- grid.arrange(gInit + marginOK, left=yTitle)
  gBEClust <- grid.arrange(gC2   + marginOK, left=yTitle)
  gCAll    <- grid.arrange(gAll  + marginOK, left=yTitle, bottom=xTitle)
  
  gFinal <- ggarrange(gAInit, gBEClust, gCAll,
                      ncol=1,nrow=3,
                      heights=c(1,1,2),
                      labels=c(paste0('All transcripts in cluster ',originalCluster),
                               'Clustering by expression',
                               'Representative individual transcripts'),
                      vjust = 1,hjust=0,
                      font.label = list(size = 8,face='bold'))

  geneExport <- data2Clust[,c('name','cluster')]
  geneExport$mainCluster <- originalCluster
  
  gRet <- list(g1 = gInit,
               g2 = gC2,
               g3 = gAll,
               gMerge = gFinal,
               geneList = geneExport,
               data = data2Clust,
               plotData = allBoth)
  
  return(gRet)
}  

for (nC in 1:optimalK){
  if (exists("c")){remove(c)}
  c <- clusterCluster(eNames = expNames,
                      sNames = stageNames,
                      data2Clust = dTSS[dTSS$newClusterNumber == nC,],
                      originalCluster = nC)

  nW <- 5 ; nH <- 7
  
  myPNG <- getIMGname(fname = paste0('Figure_S5_Cluster',nC),
                      type='PNG',
                      saveLocation = imgdir)
  
  myPDF <- getIMGname(fname = paste0('Figure_S5_Cluster',nC),
                      type='PDF',
                      saveLocation = imgdir)
  graphics.off()
  
  png(myPNG,width=nW,height=nH,res=400,units='in')
  print(ggarrange(c$gMerge))
  dev.off()
  
  pdf(myPDF,width=nW,height=nH)
  print(ggarrange(c$gMerge))
  dev.off()
  
}
