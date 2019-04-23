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
#source(paste0(workdir,'accessoryFiles/scripts/R/sortingPaper_functions.R'))
source(paste0(workdir,'accessoryFiles/scripts/R/genericFunctionsLamEtAl.R'))
source(paste0(workdir,'accessoryFiles/scripts/R/sortingPaper_plotExpressionForSingleGene_v_Bortvin.R'))

library(ggplot2)
library(reshape2)
library(plyr)

theme7point()

hstss <- read.table(paste0(outdir,
                           '/Rtables/hstssTable.forR.tab'),
                    header=TRUE)

## Hotspot, not a promoter (or putative functional site (KO hotspot))
hsH3 <- hstss[hstss$HS & !hstss$HSPrdm9KO & !hstss$TSS & hstss$spo11>0 & hstss$cs != 'chrX',]

hs1    <- hsH3[hsH3$LE.R1.H3K4me3.peaks & hsH3$ZY.R1.H3K4me3.peaks,]
hs1$LZ <- log2(hs1$LE.R1.H3K4me3.KM/hs1$ZY.R1.H3K4me3.KM)
hs1$LZ <- hs1$LZ - median(hs1$LZ,na.rm=TRUE)

hs2    <- hsH3[hsH3$LE.R2.H3K4me3.peaks & hsH3$ZY.R2.H3K4me3.peaks,]
hs2$LZ <- log2(hs2$LE.R2.H3K4me3.KM/hs2$ZY.R2.H3K4me3.KM)
hs2$LZ <- hs2$LZ - median(hs2$LZ,na.rm=TRUE)

hs3    <- hsH3[hsH3$LE.R2.abK4me.peaks & hsH3$ZY.R2.abK4me.peaks,]
hs3$LZ <- log2(hs3$LE.R2.abK4me.NCISTSS/hs3$ZY.R2.abK4me.NCISTSS)
hs3$LZ <- hs3$LZ - median(hs3$LZ,na.rm=TRUE)

nC <- 5

useOligoStrength <- TRUE

if (useOligoStrength){
  hs1$q5Strength  <- convertToQuantiles(hs1$spo11,nC)
  hs2$q5Strength  <- convertToQuantiles(hs2$spo11,nC)
  hs3$q5Strength  <- convertToQuantiles(hs3$spo11,nC)
  xLabelTxt       <- 'Hotspots (ranked by Spo11-oligo signal)'
}else{
  hs1$q5Strength  <- convertToQuantiles(hs1$HSstrength,nC)
  hs2$q5Strength  <- convertToQuantiles(hs2$HSstrength,nC)
  hs3$q5Strength  <- convertToQuantiles(hs3$HSstrength,nC)
  xLabelTxt        <- 'Hotspots (ranked by SSDS signal)'  
}

gLZ1 <- ggplot(hs1,aes(x=q5Strength,y=LZ)) + 
  geom_hline(yintercept=0,lwd=.2, color='grey40') + 
  geom_violin(fill='steelblue',
              alpha=.4,
              lwd=.2) + 
  geom_boxplot(notch=TRUE,
               width=.3,
               fill='grey99',
               alpha=.8,
               outlier.size=.1,
               lwd=.2) + 
  scale_x_discrete(breaks=c(1,nC),labels=c('Weak','Strong')) + 
  xlab(xLabelTxt) + 
  ylab('H3K4me3 ratio\nlog2 (Leptotene/Zygotene) - median')

gLZ2 <- ggplot(hs2,aes(x=q5Strength,y=LZ)) + 
  geom_hline(yintercept=0,lwd=.2, color='grey40') + 
  geom_violin(fill='steelblue',
              alpha=.4,
              lwd=.2) + 
  geom_boxplot(notch=TRUE,
               width=.3,
               fill='grey99',
               alpha=.8,
               outlier.size=.1,
               lwd=.2) + 
  scale_x_discrete(breaks=c(1,nC),labels=c('Weak','Strong')) + 
  xlab(xLabelTxt) + 
  ylab('H3K4me3 ratio\nlog2 (Leptotene/Zygotene) - median')

gLZ3 <- ggplot(hs3,aes(x=q5Strength,y=LZ)) + 
  geom_hline(yintercept=0,lwd=.2, color='grey40') + 
  geom_violin(fill='steelblue',
              alpha=.4,
              lwd=.2) + 
  geom_boxplot(notch=TRUE,
               width=.3,
               fill='grey99',
               alpha=.8,
               outlier.size=.1,
               lwd=.2) + 
  scale_x_discrete(breaks=c(1,nC),labels=c('Weak','Strong')) + 
  xlab(xLabelTxt) + 
  ylab('abK4me ratio\nlog2 (Leptotene/Zygotene) - median')

sliceFig <- drawSlice(sliceLEDZ = paste0(outdir,
                                         '/slices/LE.R1.H3K4me3_ncis.chr10_92Mb.DZ.tab'),
                      smoothingWindow = 10)
sliceFig$gFig

gBC <- ggarrange(gLZ1, gLZ3, 
                  nrow=1,ncol=2,
                  labels=c('b','c'),
                  vjust = 1,hjust=0,
                  font.label = list(size = 8,
                                    face='bold'))

gABC <- ggarrange(sliceFig$gFig,gBC, 
                 nrow=2,ncol=1,
                 labels=c('a',''),
                 heights=c(6,4),
                 vjust = 1,hjust=0,
                 font.label = list(size = 8,
                                   face='bold'))

png(getIMGname(fname = 'Figure_S3',
               type = 'PNG',
               saveLocation = imgdir),
    height=5,width=5,units='in',res=300)
print(gABC)
dev.off()

png(getIMGname(fname = 'Figure_S3',
               type = 'PDF',
               saveLocation = imgdir),
    height=5,width=5)
print(gABC)
dev.off()
