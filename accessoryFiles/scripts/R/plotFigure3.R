#Sys.setenv(KBPIPEOUTDIR='/data/RDCO/kevbrick/GL_Sorting_Paper/finalPipe/output/')
#Sys.setenv(KBPIPEWORKDIR='/data/RDCO/kevbrick/GL_Sorting_Paper/')

## Load custom functions & H3 data
workdir <- Sys.getenv('KBPIPEWORKDIR')
outdir  <- Sys.getenv('KBPIPEOUTDIR')
datadir <- paste0(outdir,'/Rtables/')

if (exists('testingKBGLPIPE')){
  imgdir = '/home/kevbrick/testImg/'
}else{
  imgdir  <- './'
}

print(workdir)

## Load custom functions & H3 data
source(paste0(workdir,'accessoryFiles/scripts/R/genericFunctionsLamEtAl.R'))
source(paste0(workdir,'accessoryFiles/scripts/R/sortingPaper_plotExpressionForSingleGene_v_Bortvin.R'))

library(ggplot2)
library(reshape2)
library(plyr)
library(scales)

theme7point()

###########################################################################
# PANEL A
#pHS  <- plotFromDeepToolsMatrix(inMatrix = paste0(outdir,'/figs/H3K4me3.R1.ncisCorrectedByStablePromoters.HS.matrix.gz'),
pHS  <- plotFromDeepToolsMatrix(inMatrix = paste0(outdir,'/figs/H3K4me3.R1.monoCorrected.HS.matrix.gz'),
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
  ylab('H3K4me3 signal (A.U.)') + facet_grid(~stage)

myHM <- pHS$figHM + 
  scale_fill_gradientn(colors=c('white','grey80','black','red','firebrick'),values=c(0,.2,.5,.7,1)) + 
  ylab('H3K4me3 signal\nat DSB hotspots') + 
  scale_x_continuous(breaks=c(-2,0,2)) + 
  coord_cartesian(xlim=c(-2.5,2.5))

gA <- ggarrange(pHS$figMean, myHM + theme(strip.text=element_blank()), gBoxPlot ,
                ncol=1,nrow=3,
                heights=c(2,5,5),
                align='v',
                font.label = list(size = 8,face='bold'),
                vjust=1,hjust=0,
                labels=c('a',''))

###########################################################################
# PANEL B; CoveragePlot
sliceFigs <- drawSlice(sliceLEDZ = paste0(outdir,
                                          '/slices/LE.R1.H3K4me3_ncis.chr1_177Mb.DZ.tab'),
                                          smoothingWindow=10, winSz=1, RPK=FALSE)
gB <- sliceFigs$gFig

###########################################################################
# PANEL C; H3K4me3 v Expression
###########################################################################
source(paste0(workdir,'accessoryFiles/scripts/R/plotFigure3_Final_H3K4me3VExpression.R'))
gC <- gK4vExp

###########################################################################
# PANEL ABCD; 

gAB  <- ggarrange(gA,gB,
                   ncol=2,nrow=1,
                   widths=c(3,4),
                   font.label = list(size = 8,face='bold'),
                   vjust=1,hjust=0,
                   labels=c('','b'))

gABC <- ggarrange(gAB,gC,
                  ncol=1,nrow=2,
                  heights=c(5,3),
                  font.label = list(size = 8,face='bold'),
                  vjust=1,hjust=0,
                  labels=c('','c'))

ggsave(plot=gABC, 
       filename = getIMGname(saveLocation = imgdir,
                             fname = 'Figure3',
                             type = 'PNG'), 
       width=5, height=5)

ggsave(plot=gABC, 
       filename = getIMGname(saveLocation = imgdir,
                             fname = 'Figure3',
                             type = 'PDF'), 
       width=5, height=5)

###########################################################################
# Plot supplementary figure(s) related to figure 3; 

## Fig. S2
source(paste0(workdir,'accessoryFiles/scripts/R/sortingPaper_plotSupplementaryFigure3_H3K4me3_at_strongestHS.R'))

## Fig. S3
source(paste0(workdir,'accessoryFiles/scripts/R/sortingPaper_clustering_of_clusters_forSuppFig.R'))

## Fig. S4
source(paste0(workdir,'accessoryFiles/scripts/R/sortingPaper_plotUnbiasedClusteringSupplement.R'))

## Fig. S4abK4me
source(paste0(workdir,'accessoryFiles/scripts/R/sortingPaper_plotUnbiasedClusteringSupplementabK4me.R'))
