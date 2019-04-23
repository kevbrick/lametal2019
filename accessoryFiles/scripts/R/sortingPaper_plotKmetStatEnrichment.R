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
library(data.table)

theme7point()

#########################################################
df           <- fread('KmetStat_Enrichment.tab',header=TRUE)

totDF        <- aggregate(df$n,
                          by=list(stage=df$stage,rep=df$rep,ab=df$ab),
                          FUN=sum)

dfT          <- plyr:::join(df,totDF,by=c('stage','rep','ab'))

ddf          <- plyr:::join(dfT,dfT[dfT$ab=='input',],by=c('stage','rep','si'))

names(ddf)   <- c(names(ddf)[1:6],"i","nI","totI")

ddf    <- ddf[grep("_A",ddf$si),]
ddf$si <- sub("_[AB]","",ddf$si)

ddf <- ddf[ddf$ab != 'H3K9Ac',]

ddf$chipFrac <- ddf$n/ddf$x
ddf$inFrac   <- ddf$nI/ddf$totI

ddfOK        <- ddf[ddf$ab != 'input' & !(ddf$ab=='H3K9Ac' & ddf$rep == 'R2'),]
ddfOK$ratio  <- ddfOK$chipFrac/ddfOK$inFrac

dAgg         <- aggregate(ddfOK$ratio,by=list(ddfOK$stage, ddfOK$rep,ddfOK$ab,ddfOK$si),FUN=mean)

names(dAgg) <- c('stage','rep','ab','si','ratio')

dAgg$stage <- factor(dAgg$stage,levels=c('LE','ZY','EP','LP','DI'))
dAgg$ab    <- factor(dAgg$ab,levels=c('H3K4me3','abK4me'))

dAgg$lbl <- ''
dAgg$lbl[dAgg$ratio>1] <- paste0(' ',round(dAgg$ratio[dAgg$ratio>1],0))

dAgg$si[dAgg$si == 'WT'] <- 'Unmodified'

dAgg$sName <- paste0(dAgg$ab," ( ",dAgg$rep," )")
dAgg$sName[dAgg$sName == 'abK4me ( R2 )'] <- 'abK4me'
dAgg$sName    <- factor(dAgg$sName,levels=c('H3K4me3 ( R1 )','H3K4me3 ( R2 )','abK4me'))

ggplot(dAgg,aes(x=si,y=ratio,color=ab)) + 
  geom_bar(stat='identity',width=.5) +
  geom_text(aes(label=lbl,hjust=0),
            size=7*5/14)+
  coord_flip(ylim = c(0,20)) + 
  facet_grid(stage~sName) + 
  theme(legend.position='none',
        panel.spacing.x = unit(0,'cm'),
        panel.grid.major.y = element_line(size=0.2,color='grey80'),
        strip.background=element_blank(),
        strip.text = element_text(face = 'bold',size=8)) + 
  ylab('Spike-in enrichment (fold over input)') + xlab('Nucleosome type') + 
  scale_color_manual(values=c('darkorange','steelblue'))

ggsave(getIMGname(fname = 'Figure_SX_KmetStat',
                  type = 'PNG',
                  saveLocation = imgdir),
       height=8,width=6)

ggsave(getIMGname(fname = 'Figure_SX_KmetStat',
                  type = 'PDF',
                  saveLocation = imgdir),
       height=8,width=6)