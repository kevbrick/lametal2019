#Sys.setenv(KBPIPEOUTDIR='/data/RDCO/RDCO/kevbrick/GL_Sorting_Paper/forPaper/output/')
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
#source(paste0(workdir,'accessoryFiles/scripts/R/sortingPaper_functions.R'))
source(paste0(workdir,'accessoryFiles/scripts/R/genericFunctionsLamEtAl.R'))
source(paste0(workdir,'accessoryFiles/scripts/R/sortingPaper_plotExpressionForSingleGene_v_Bortvin.R'))

library(ggplot2)
library(reshape2)
library(plyr)
library(data.table)
library(ggpubr)

theme7point()

pH3R1   <- plotFromDeepToolsMatrix(inMatrix = paste0('/data/RDCO/kevbrick/GL_Sorting_Paper/finalPipe/output/figs/H3K4me3.R1.monoCorrected.HS.matrix.gz'),                  withHeatmap = TRUE, matrixType = 'referencepoint',normalizeByIgG = FALSE, sortOrder = c('LE','ZY','EP','LP','DI'))
pH3R1p  <- plotFromDeepToolsMatrix(inMatrix = paste0('/data/RDCO/kevbrick/GL_Sorting_Paper/finalPipe/output/figs/H3K4me3.R1.ncisCorrectedByStablePromoters.HS.matrix.gz'), withHeatmap = TRUE, matrixType = 'referencepoint',normalizeByIgG = FALSE, sortOrder = c('LE','ZY','EP','LP','DI'))
pH3R2   <- plotFromDeepToolsMatrix(inMatrix = paste0('/data/RDCO/kevbrick/GL_Sorting_Paper/finalPipe/output/figs/H3K4me3.R1.monoCorrected.HS.matrix.gz'),                  withHeatmap = TRUE, matrixType = 'referencepoint',normalizeByIgG = FALSE, sortOrder = c('LE','ZY','EP','LP','DI'))
pH3R2p  <- plotFromDeepToolsMatrix(inMatrix = paste0('/data/RDCO/kevbrick/GL_Sorting_Paper/finalPipe/output/figs/H3K4me3.R1.ncisCorrectedByStablePromoters.HS.matrix.gz'), withHeatmap = TRUE, matrixType = 'referencepoint',normalizeByIgG = FALSE, sortOrder = c('LE','ZY','EP','LP','DI'))
pAb     <- plotFromDeepToolsMatrix(inMatrix = paste0('/data/RDCO/kevbrick/GL_Sorting_Paper/finalPipe/output/figs/abK4me.R2.ncisCorrectedByStablePromoters.HS.matrix.gz'),  withHeatmap = TRUE, matrixType = 'referencepoint',normalizeByIgG = FALSE, sortOrder = c('LE','ZY','EP','LP','DI'))

getFig <- function(x,xTit){
  return(ggarrange(x$figMean,
                   x$figHM + theme(strip.text = element_blank()) + 
                     ylab('Hotspots'),
                   ncol=1,nrow=2,
                   heights=c(4,9),
                   align='v') + ggtitle(xTit) + theme(plot.title = element_text(hjust=.5,size=8)))
}

gH3R1  <- getFig(pH3R1,"H3K4me3 (R1; KM)")
gH3R1p <- getFig(pH3R1p,"H3K4me3 (R1; TSS)")
gH3R2  <- getFig(pH3R2,"H3K4me3 (R2; KM)")
gH3R2p <- getFig(pH3R2p,"H3K4me3 (R2; TSS)")
gAb    <- getFig(pAb,"H3K4me2/3 (TSS)")

png(getIMGname('Figure_S4','PNG',
               saveLocation = imgdir),
    width=7,height=9,
    units='in',res=300);
print(ggarrange(gH3R1p, 
                gH3R1, 
                gH3R2p,
                gH3R2, 
                gAb, 
                labels=c('a','b','c','d','e',''),
                vjust = 1,hjust=0,
                font.label = list(size = 8,
                                  face='bold'),
                ncol=2,nrow=3))
dev.off()
