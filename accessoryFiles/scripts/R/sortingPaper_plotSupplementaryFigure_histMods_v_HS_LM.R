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

library(gridExtra)
library(grid)
library(ggfortify)
library(ggplot2)
library(factoextra)
library(leaps)
library(numform)
library(ggpubr)
library(corrplot)
library(png)

## Functions
plotMLR <- function(r=NULL,
                    nProp=7,
                    name='generic') {
  
  ## get Best N subsets
  leaps <- regsubsets(as.formula(paste("strength ~ ", 
                                       paste(names(r)[1:nProp], 
                                             collapse=" + "))), 
                      data = r,
                      nbest=1,
                      nvmax=nProp)
  
  # Get Summary
  s <- summary(leaps)
  
  # make R2 v N DF
  dR2 <- data.frame(R2=s$rsq,
                    diff=s$rsq-c(0,s$rsq[2:nProp]),
                    N=1:nProp)
  
  # Melt to usable DF for ggplot
  dfS <- cbind(data.frame(N=1:nProp),
               as.data.frame(s$which))
  
  # Get the best order (most NB first)
  m <- reshape2:::melt.data.frame(dfS,id.vars = "N")
  
  m <- m[m$variable != '(Intercept)',]
  dd <- as.data.frame(table(m[,2:3]))
  dd <- dd[dd$value==TRUE,]
  dd <- dd[order(dd$Freq),]
  
  m$variable <- factor(m$variable,
                       levels = dd$variable)
  
  m$lbl          <- ''
  m$lbl[m$value] <- 'x'
  
  yMin <- min(dR2$R2)
  yMax <- max(dR2$R2)
  
  yMin <- yMin - (yMax-yMin)*0.1
  yMax <- yMax + (yMax-yMin)*0.1
  
  gR2 <- ggplot(dR2,aes(x=N,y=R2)) + 
    geom_line(lwd=0.2) + 
    geom_point(size=0.4,lwd=.2,shape=21,fill='orange') + 
    xlab('Histone marks in model') + 
    ylab(expression('Correlation with SSDS (R'^2*')')) +
    theme()  +
    scale_x_continuous(breaks=c(1,3,5,7)) +
    coord_cartesian(xlim=c(0.5,nProp+0.5),
                    ylim=c(yMin,yMax),
                    expand = FALSE)
  
  gDiff <- ggplot(dR2,aes(x=N,y=diff)) + 
    geom_line(lwd=.2) + 
    geom_point(size=0.4,lwd=.2,shape=21,fill='green') + 
    xlab('') + 
    ylab(expression('R'^2*' improvement')) +
    theme(axis.text.x=element_blank()) +
    coord_cartesian(xlim=c(0.5,nProp+0.5),expand = FALSE)
  
  gProp <- ggplot(m,aes(y=variable,x=N,fill=value)) + 
    geom_tile(size=0.2,color='black') + 
    scale_fill_manual(values=c("white","salmon")) + 
    theme(panel.grid = element_blank(),
          legend.position='none') + 
    ylab('') + 
    xlab('Model features (#)') + 
    geom_text(aes(label=lbl),
              size=7*5/14) +
    theme(axis.ticks.y=element_blank(),
          axis.ticks.x=element_blank(),
          axis.line = element_blank()) + 
    scale_x_continuous(breaks=c(1,3,5,7)) + 
    coord_cartesian(xlim=c(0.5,nProp+0.5),expand = FALSE)
  
  gBoth <- ggarrange(gR2 + theme(axis.text.x=element_blank()) + xlab(''),gProp,
                     ncol=1,nrow=2,
                     heights=c(2,1),
                     align='v')
  
  return(list(ggR2=gR2,
              ggRUp=gDiff,
              gGProp=gProp,
              ggAll=gBoth))
}

hstss <- read.table(paste0(outdir,
                           '/Rtables/hstssTable.forR.tab'),
                    header=TRUE)

## remove duplicates
## Pick subset of hist. mods. 
col2Use <- names(hstss) %in% c('H3K4me3','H3K9ac','H3K4me2','H3K4me1','H3K36me3','H3K27ac','H4ac5','HSstrength')
row2Use <- hstss$cs %in% paste0("chr",1:19) & hstss$HS & !hstss$TSS & !hstss$HSPrdm9KO

iData <- hstss[row2Use,col2Use] 
names(iData)[names(iData) == 'HSstrength'] <- 'strength'

columnOrder <- c('H3K4me3','H3K9ac','H3K4me2','H3K4me1','H3K36me3','H3K27ac','H4ac5','strength');
iData <- iData %>% select(.,columnOrder)

iData$noZeros <- apply(iData[,1:8],1,FUN=function(x){sum(x<=0)}) < 1
  
iOK <- iData[iData$noZeros,1:8]

## Plot CC Matrix
ccPlot <- ggCorMat(cor(iOK)^2,
                   decimalPlaces = 2,
                   keepLeadingZeros = TRUE,
                   flipIt = TRUE) + 
  scale_y_discrete(position = 'right') + 
  geom_vline(xintercept=7.5)

## Plot LM Figure
ggMLR <- plotMLR(iOK,7)

## Join figs & export
ggCorrelationPlot <- ggarrange(ccPlot,ggMLR$ggAll,
                               ncol=2,nrow=1,
                               labels=c('a','b'),
                               vjust = 1,hjust=0,
                               font.label = list(size = 8,face='bold')
)

## Plot fig S10
png(file = getIMGname(saveLocation = imgdir,
                      fname = 'Figure_S10',
                      type = 'PNG'),
    width=7,height=3,
    units='in',res=300)
print(ggCorrelationPlot)
dev.off()

pdf(file = getIMGname(saveLocation = imgdir,
                      fname = 'Figure_S10',
                      type = 'PDF'),
    width=7,height=3)
print(ggCorrelationPlot)
dev.off()