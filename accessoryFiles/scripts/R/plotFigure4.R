## This script plots figure 4 and Figs. S9-11

#Sys.setenv(KBPIPEOUTDIR='/data/RDCO/kevbrick/GL_Sorting_Paper/finalPipe/output/')
#Sys.setenv(KBPIPEWORKDIR='/data/RDCO/kevbrick/GL_Sorting_Paper/')

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

library(png)
library(ggpubr)
library(plyr)

## Functions
plot3A <- function(sNum=1,
                   featureSize=1000,
                   xlimits=NULL){
  grHS  <- read.table(paste0(outdir,'/slices/fig4/slice',sNum,'.HSB6.bed'))
  grHS  <- grHS[,1:3]
  names(grHS) <- c('cs','from','to')
  
  grHS$mid  <- round((grHS$from+grHS$to)/2)
  
  if (is.null(featureSize)){
    grHS$minX <- grHS$from
    grHS$maxX <- grHS$to
  }else{
    grHS$minX <- grHS$mid - featureSize
    grHS$maxX <- grHS$mid + featureSize
  }
  
  grHS$minY <- 1
  grHS$maxY <- 2
  grHS$type <- 'Hotspots'
  
  ################################
  grTSS  <- read.table(paste0(outdir,'/slices/fig4/slice',sNum,'.TSSgencode.bed'))
  grTSS  <- grTSS[,1:3]
  names(grTSS) <- c('cs','from','to')
  grTSS$mid  <- round((grTSS$from+grTSS$to)/2)
  
  if (is.null(featureSize)){
    grTSS$minX <- grTSS$from
    grTSS$maxX <- grTSS$to
  }else{
    grTSS$minX <- grTSS$mid - featureSize
    grTSS$maxX <- grTSS$mid + featureSize
  }
  
  grTSS$minY <- 1
  grTSS$maxY <- 2
  grTSS$type <- 'TSS'
  
  cAnnotation <- rbind(grHS,grTSS)
  
  ################################
  cAll <- read.table(paste0(outdir,'/slices/fig4/slice',sNum,'.ALL.bedgraph'))
  
  names(cAll) <- c('cs','pos','cover','type')
 
  ## Do not use H3K4me3 from rebuilt pool
  cAll <- cAll[cAll$type != 'abK4me',]
  
  cAgg <- aggregate(cAll$cover, by=list(cAll$type),FUN=sum, na.rm=TRUE)
  names(cAgg) <- c('type','tot')
  
  cBoth <- join(x=cAll,y=cAgg,by='type')
  
  cBoth$TPM <- cBoth$cover/cBoth$tot*1000000
  
  cAll <- cBoth
  
  cAll$type <- factor(cAll$type,levels=c('H3K4me3',
                                         'H3K9ac',
                                         'H3K36me3',
                                         'H3K4me2',
                                         'H3K4me1',
                                         'H3K27ac',
                                         'H4ac5',
                                         'H4K20me3',
                                         'H4K12ac',
                                         'H4K8ac',
                                         'H3K79me1',
                                         'H3K27me3',
                                         'H3K79me3',
                                         'H3K4ac',
                                         'H3K27me1',
                                         'H3K9me3',
                                         'H3K9me2',
                                         'H3',
                                         'IgG',
                                         'Input'))
  
  getMaxYLim <- function(xx){
    maxVals    <- c(50,100,500,1000,3000,10000,100000)
    mV         <- (maxVals-xx)
    mV[mV < 0] <- 999999
    myMax      <- maxVals[which(mV==min(mV))[1]]
    return(myMax)
  }
  
  if (is.null(xlimits)){
    xAllMin <- min(cAll$pos)/1000000
    xAllMax <- max(cAll$pos)/1000000
  }else{
    xAllMin <- xlimits[1]
    xAllMax <- xlimits[2]
  }
  
  aggDet <- aggregate(cAll$TPM,
                  by=list(cAll$type),
                  FUN=function(x){round(max(x))})
  names(aggDet) <- c('type','maxVal')
  
  aggDet$n <- 1:length(aggDet$max)
  
  aggDet$maxComma <- format(aggDet$max, big.mark=",")
  aggDet$figName  <- paste0(aggDet$type,' : ',aggDet$maxComma)
  
  cAll <- join(cAll,aggDet,by='type')
  
  cAll$figName <- factor(cAll$figName,
                         levels=aggDet$figName)
  
  gInit <- ggplot(cAll,
                  aes(x=pos/1000000,
                      y=TPM,
                      fill=type,
                      color=type)) 
  
  
  for (hs in 1:dim(grHS)[1]){
    gHSfill <-annotate(geom='rect',
                       xmin=grHS$minX/1000000,
                       xmax=grHS$maxX/1000000,
                       ymin=-Inf,
                       ymax=Inf,
                       fill=alpha('forestgreen',.075),
                       color=NA)
  }
  
  for (hs in 1:dim(grTSS)[1]){
    gTSSfill <- annotate(geom='rect',
                         xmin=grTSS$minX/1000000,
                         xmax=grTSS$maxX/1000000,
                         ymin=-Inf,
                         ymax=Inf,
                         fill=alpha('orange',.075),
                         color=NA)
  }
  
  #scale_color_manual(values=getHistoneModColorsColors()) + 
  #scale_fill_manual(values=getHistoneModColorsColors()) + 
  g3ACoverage <- gInit + 
    geom_polygon(lwd=0.1,fill='grey50',color='black') + 
    facet_wrap(~type,ncol=1,strip.position = 'left') + 
    gHSfill + gTSSfill +
    xlab(paste0('Position on ',cAll$cs[1],' (Mbp)')) + 
    ylab('Coverage') + 
    theme(legend.position='none',
          strip.background = element_blank(),
          strip.text.y = element_text(size=7,angle = 180),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank(),
          axis.line.y=element_blank(),
          axis.line.x=element_blank(),
          plot.margin=unit(c(0,0,0,0),'cm'),
          panel.spacing = unit(0,'cm')) + 
    coord_cartesian(xlim=c(xAllMin,xAllMax),expand = FALSE)
  
  g3ACoverageFreeY <- gInit + 
    geom_polygon(lwd=0.3,
                 fill='grey50',color='black') + 
    facet_grid(figName~.,
               scales='free_y') +
    gHSfill + gTSSfill +
    xlab(paste0('Position on ',cAll$cs[1],' (Mbp)')) + 
    ylab('Coverage') + 
    theme(legend.position='none',
          strip.background = element_blank(),
          strip.text.y = element_text(angle = 0,hjust = 1),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank(),
          axis.line.y=element_blank(),
          axis.line.x=element_blank(),
          plot.margin=unit(c(0,0,0,0),'cm')) + 
    coord_cartesian(xlim=c(xAllMin,xAllMax),expand = FALSE)

  gAnnotation <- ggplot(cAnnotation,
                        aes(xmin=minX/1000000,
                            xmax=maxX/1000000,
                            ymin=minY,
                            ymax=maxY,
                            fill=type,
                            color=type)) + 
    geom_rect(lwd=0.2) + 
    facet_wrap(~type,nrow=2,ncol=1,strip.position='left') + 
    coord_cartesian(xlim=c(xAllMin,xAllMax),ylim=c(0.95,2.05),expand = FALSE) +
    scale_color_manual(values=c(alpha(darkenColor('forestgreen',1),1),
                                alpha(darkenColor('orange',1))),1) + 
    scale_fill_manual(values=c(alpha(darkenColor('forestgreen',1.3),1),
                               alpha(darkenColor('orange',1.3),1))) + 
    xlab(paste0('Position on ',cAll$cs[1],' (Mbp)')) + 
    ylab('') + 
    theme(legend.position='none',
          axis.text.y  = element_blank(),        
          axis.ticks.y = element_blank(),
          axis.line.y = element_blank(),
          strip.background = element_blank(),
          strip.text.y = element_text(size=7,angle=180),
          plot.margin = unit(c(-1,0,0,0),'cm')) 
    
  # geom_text(aes(x, y, label=lab),
  #             data=data.frame(x=Inf, y=1.5, 
  #                             minX=Inf,maxX=Inf,
  #                             minY=Inf,maxY=Inf,
  #                             lab=c("Hotspots","TSS"),
  #                             type=c("Hotspots","TSS")),
  #             vjust=0.5, hjust=1, size=7*5/14)
  
  themeNoX    <- theme(axis.text.x = element_blank())
  themeNoMarg <- theme(plot.margin=unit(c(0,0,0,0),'cm'))
  
  gCoverage <- ggarrange(g3ACoverage + xlab('') + themeNoX + themeNoMarg,
                         gAnnotation + themeNoMarg, 
                         align = 'v',
                         heights=c(14,2),
                         ncol=1,nrow=2)
  
  gCoverageFree <- ggarrange(g3ACoverageFreeY + themeNoX + themeNoMarg,
                             gAnnotation + themeNoMarg, 
                         align = 'v',
                         heights=c(14,2),
                         ncol=1,nrow=2)
  
  return(list(figRPKM=gCoverage,
              figFree=gCoverageFree))
  
}

## Panel A
#gAx2 <- plot3A(sNum=2)
gAx2 <- plot3A(sNum=4, xlimits=c(31.45,31.525))
g4A <- gAx2$figRPKM


## Panel B
# NOTE: Also plots Supp. Fig 9
source(paste0(workdir,'accessoryFiles/scripts/R/sortingPaper_histoneModHeatMaps_atGenomicFeatures.R'))
g4B     <-pScaleHS$fig + 
  geom_vline(xintercept=50,color='grey50',alpha=.3) + 
  scale_fill_gradient2(low='royalblue1',mid='white',high='red',midpoint = 0)+
  scale_x_continuous(breaks=c(0,50,100),labels=c(-5,0,5)) + 
  coord_cartesian(xlim=c(0,100),expand = TRUE) + 
  theme(legend.title=element_text(),
        legend.key.width = unit(0.25,'cm'),
        legend.spacing = unit(c(0,0,0,0),'cm'),
        legend.box.margin = unit(c(0,0,0,-0.2),'cm'),
        legend.direction='vertical') + 
  xlab('Distance to hotspot center (Kbp)')

###########################################################################
# PANEL C
pHS  <- plotFromDeepToolsMatrix(inMatrix = paste0(outdir,'/figs/H3K9Ac.R1.ncisCorrectedByStablePromoters.HS.matrix.gz'),
                                matrixType = 'referencepoint',
                                normalizeByIgG = FALSE,
                                normByRow = FALSE,
                                withHeatmap = TRUE,
                                sortOrder = c('LE','ZY','EP','LP','DI'))

gBoxPlot <- ggplot(pHS$dDist,aes(x=1,y=mid/flank,fill=stage)) + 
  geom_hline(lwd=.3,lty='dashed',color='magenta',yintercept=.75) + 
  geom_violin(lwd=.1) + scale_y_log10(label=fancy_scientific) + 
  geom_boxplot(width=.3,fill='grey90',color='black',lwd=.1,
               outlier.size=.1, 
               outlier.alpha=.1,
               notch=TRUE) + 
  coord_cartesian(ylim=c(1e-3,1e4)) + 
  scale_fill_manual(values=getStageColors()) + 
  scale_color_manual(values=getStageColors(darkLvl = 1.8)) + 
  theme(legend.position='none',
        strip.background=element_blank(),
        strip.text=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x = element_blank()) + 
  xlab('Meiosis stage') + 
  ylab('H3K9ac signal (A.U.)') + facet_grid(~stage)

myHM <- pHS$figHM + 
  scale_fill_gradientn(colors=c('white','grey80','black','red','firebrick'),values=c(0,.2,.5,.7,1)) + 
  ylab('H3K9ac signal\nat DSB hotspots') + 
  scale_x_continuous(breaks=c(-2,0,2)) + 
  coord_cartesian(xlim=c(-2.5,2.5))

g4C <- ggarrange(pHS$figMean, myHM + theme(strip.text=element_blank()), gBoxPlot ,
                ncol=1,nrow=3,
                heights=c(2,5,5),
                align='v',
                font.label = list(size = 8,face='bold'),
                vjust=1,hjust=0,
                labels=c('c',''))

## PANEL D
## NOTE: Also plot figS11
source(paste0(workdir,'accessoryFiles/scripts/R/sortingPaper_histoneModPCA_vFinal.R'))

g4D <- gROCs$gAUCs

g4E <- gContribs

g4B <- g4B + theme(legend.position='top',
                   legend.direction = 'horizontal',
                   legend.key.width=unit(0.6,'cm'),
                   legend.key.height=unit(0.2,'cm'))

g4ABC <- ggarrange(g4A,g4B,g4C,ncol=3,nrow=1,
                   labels=c('a','b','c'),
                   vjust = 1,hjust=0,
                   font.label = list(size = 8,face='bold')) 

g4DE <- ggarrange(g4D,g4E,ncol=2,nrow=1,
                   labels=c('d','e'),
                   vjust = 1,hjust=0,
                   font.label = list(size = 8,face='bold')) 


g4All <- ggarrange(g4ABC,g4DE,ncol=1,nrow=2,
                   heights=c(6,4),
                   vjust = 1,hjust=0,
                   font.label = list(size = 8,face='bold'))


png(file = getIMGname(saveLocation = imgdir,
                      fname = 'Figure4',
                      type = 'PNG'),
    width=7,height=6,
    units='in',res=300)
print(g4All)
dev.off()

pdf(file = getIMGname(saveLocation = imgdir,
                      fname = 'Figure4',
                      type = 'PDF'),
    width=7,height=6)
print(g4All)
dev.off()

#####################
## Plot S10
source(paste0(workdir,'accessoryFiles/scripts/R/sortingPaper_plotSupplementaryFigure_histMods_v_HS_LM.R'))