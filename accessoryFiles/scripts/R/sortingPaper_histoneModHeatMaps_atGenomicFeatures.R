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
  imgdir <- './'
}

## Load custom functions & H3 data
source(paste0(workdir,'accessoryFiles/scripts/R/genericFunctionsLamEtAl.R'))

library(ggplot2)
library(reshape2)
library(plyr)

theme7point()

normByIgG <- TRUE
theme7point() ## SET THEME

SO <- c("abK4me","H3K4me3" ,"H3K9ac"  ,"H3K36me3","H3K4me2","H3K4me1",
       "H3K27ac" ,"H4ac5"   ,"H4K20me3","H4K12ac","H4K8ac",
       "H3K79me1","H3K27me3","H3K79me3","H3K4ac" ,"H3K27me1",
       "H3K9me3" ,"H3K9me2" ,"H3"      ,"Input"  ,"IgG")

pScaleHS  <- plotFromDeepToolsMatrix(inMatrix = paste0(outdir,'/tables/hmHS.matrix.gz'),
                                     matrixType = 'referencepoint',
                                     normalizeByIgG = normByIgG,
                                     sortOrder = SO,
                                     samplesToExclude = 'abK4me')

pScaleTSS  <- plotFromDeepToolsMatrix(inMatrix = paste0(outdir,'/tables/hmTSS.matrix.gz'),
                                     matrixType = 'referencepoint',
                                     normalizeByIgG = normByIgG,
                                     sortOrder = SO,
                                     samplesToExclude = 'abK4me')

pScaleTES  <- plotFromDeepToolsMatrix(inMatrix = paste0(outdir,'/tables/hmTES.matrix.gz'),
                                      matrixType = 'referencepoint',
                                      normalizeByIgG = normByIgG,
                                      sortOrder = SO,
                                      samplesToExclude = 'abK4me')

pScaleEnhancer  <- plotFromDeepToolsMatrix(inMatrix = paste0(outdir,'/tables/hmENH.matrix.gz'),
                                      matrixType = 'referencepoint',
                                      normalizeByIgG = normByIgG,
                                      sortOrder = SO,
                                      samplesToExclude = 'abK4me')

pScaleGene <- plotFromDeepToolsMatrix(inMatrix=paste0(outdir,'/tables/hmCDS.matrix.gz'),
                                        matrixType = 'scaleregions',
                                        normalizeByIgG = normByIgG,
                                        sortOrder=SO,
                                        samplesToExclude = 'abK4me')
  
pScaleLINE <- plotFromDeepToolsMatrix(inMatrix=paste0(outdir,'/tables/hmLINE.matrix.gz'),
                                      matrixType = 'scaleregions',
                                      normalizeByIgG = normByIgG,
                                      sortOrder=SO,
                                      samplesToExclude = 'abK4me')

noLeg <- theme(legend.position='top', 
               legend.title=element_blank(), 
               plot.title=element_text(hjust = 0.5,
                                       lineheight=.8,
                                       size=8),
               legend.key.height = unit(0.4,'cm'),
               legend.key.width = unit(0.3,'cm'),
               legend.margin=margin(c(0,0,0,0)))
noTxt <- theme(axis.text.y = element_blank())

gAll <- ggarrange(pScaleHS$fig +noLeg + ggtitle('HS'),
                  pScaleTSS$fig+noLeg+noTxt + ggtitle('TSS'),
                  pScaleTES$fig+noLeg+noTxt + ggtitle('TES'),
                  pScaleGene$fig+noLeg+noTxt + ggtitle('Gene'),
                  pScaleLINE$fig+noLeg+noTxt + ggtitle('LINEs'),
                  pScaleEnhancer$fig+noLeg+noTxt + ggtitle('Enhancer'),
                  ncol=6,nrow=1,
                  widths=c(1.4,1,1,1,1,1,1))


png(file = getIMGname(saveLocation = imgdir,
                      fname = 'Figure_S9',
                      type = 'PNG'),
    width=7,height=3,
    units='in',res=300)
print(gAll)
dev.off()

pdf(file = getIMGname(saveLocation = imgdir,
                      fname = 'Figure_S9',
                      type = 'PDF'),
    width=7,height=3)
print(gAll)
dev.off()
