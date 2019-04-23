Sys.setenv(KBPIPEOUTDIR='/data/RDCO/kevbrick/GL_Sorting_Paper/forPaper/output/')
Sys.setenv(KBPIPEWORKDIR='/data/RDCO/kevbrick/GL_Sorting_Paper/')

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

library(gridExtra)
library(grid)
library(ggfortify)
library(ggplot2)
library(ggcorrplot)
library(ggpubr)
library(factoextra)
library(leaps)
library(numform)
library(corrplot)
library(ggrepel)

## makeROC function
makeROCsPlot <- function(histmods,dataHM,dataPCA){
  
  nHM <- 0
  
  for (hm in histmods){
    nHM = nHM+1
    print(hm)
    roc <- makeROC(score = dataHM[[hm]],
                   type  = dataHM$name %in% c('HS'),
                   tryFR = TRUE )
    
    roc$data <- roc$data[seq(1,dim(roc$data)[1],20),]
    
    roc$data$type  <- hm
    roc$data$type2 <- 'HM'
    
    myAUC <- data.frame(type   = 'HistoneMod',
                        detail = hm,
                        AUC    = roc$manualAUC)
    
    if (nHM == 1){
      rocAll <- roc$data
      retAUC <- myAUC
    }else{
      rX <- rocAll
      rocAll <- rbind(rX,roc$data)
      
      rA <- retAUC
      retAUC <- rbind (rA,myAUC)
    }
  }
  
  for (pc in 1:(dim(dataPCA)[2]-1)){
    print(paste0('PC',pc))
    roc <- makeROC(score = dataPCA[[paste0('Dim.',pc)]],
                   type  = dataPCA$name %in% c('HS'),
                   tryFR = TRUE )
    
    roc$data <- roc$data[seq(1,dim(roc$data)[1],20),]
    
    roc$data$type  <- paste0('PC',pc)
    roc$data$type2 <- 'PC'
    myAUC <- data.frame(type   = 'PC',
                        detail = paste0('PC ',pc),
                        AUC    = roc$manualAUC)
    rX <- rocAll
    rocAll <- rbind(rX,roc$data)
    
    rA <- retAUC
    retAUC <- rbind (rA,myAUC)
  }
  
  rocDef <- data.frame(score=c(0,0),
                       type=c('mid','mid'),
                       TPR=c(0,1),
                       FPR=c(0,1),
                       type2=c('null','null'))
  
  rocAll$newName <- as.character(rocAll$type2)
  rocAll$newName[rocAll$type == 'PC1'] <- 'PC1'
  rocAll$newName[rocAll$type == 'H3K36me3'] <- 'H3K36me3'
  rocAll$newName[rocAll$type == 'H3K4me1'] <- 'H3K4me1'
  rocAll$newName[rocAll$type == 'H3K4me2'] <- 'H3K4me2'
  rocAll$newName[rocAll$type == 'H3K4me3'] <- 'H3K4me3'
  rocAll$newName[rocAll$type == 'H3K9ac'] <- 'H3K9ac'
  
  rocAll$type <- factor(rocAll$type,levels=c(unique(rocAll$type)[1:17],paste0('PC',1:15)))
  
  rocOK <- rocAll[rocAll$type %in% c('PC1','H3K4me3','H3K4me2','H3K36me3','H3K9ac'),]
  
  dfLabels <- data.frame(FPR=c(0.03, 0.09, 0.17, 0.5, 0.45),
                         TPR=c(1.02, 0.87, 0.5 , 0.55, 0.75),
                         lbl=c('H3K36me3','H3K4me2','PC1','H3K9ac','H3K4me3'),
                         type2=c('HM','HM','PC','HM','HM'))
  
  gNewAUCs <- ggplot(rocAll,aes(x=FPR*100,y=TPR*100,
                                group=type,color=type2)) + 
    geom_line(alpha=.2) + 
    geom_line(data=rocOK,aes(linetype=type)) + 
    geom_line(data=rocDef,color='black') + 
    geom_label(data=dfLabels,
                     aes(label=lbl,group=lbl),
                     size=7*5/14,
                     hjust=-.2) + 
    ylab('Sensitivity (hotspots; %)') + 
    xlab('1 - Specificity (hotspots; %)') + 
    theme(legend.title=element_blank(),
          legend.position=c(1,0),
          legend.justification=c(1,0),
          legend.background = element_rect(fill=NA),
          legend.spacing.y = unit('0.015','cm'),
          legend.key.height = unit('0.175','cm')) + 
    scale_color_manual(values=c('grey50','royalblue2')) 
  
  gIndividualAUCs <- ggplot(rocAll,aes(x=FPR*100,y=TPR*100,
                                group=type,color=type2)) + 
    annotate(geom='line', x=c(0,100),y=c(0,100),lty='dashed',lwd=.2) + 
    geom_line(alpha=1,lwd=.4) + 
    ylab('Sensitivity (hotspots; %)') + 
    xlab('1 - Specificity (hotspots; %)') + 
    theme(legend.title=element_blank(),
          legend.position="none",
          legend.justification=c(1,0),
          legend.background = element_rect(fill=NA),
          legend.spacing.y = unit('0.015','cm'),
          legend.key.height = unit('0.175','cm'),
          strip.text=element_blank(),
          strip.background=element_blank()) + 
    scale_color_manual(values=c('grey50','royalblue2')) +
    facet_wrap(~type) + 
    geom_text(aes(label=paste0(type,"  "),x=100,y=10,hjust=1),
              check_overlap = TRUE,
              size=7*5/14) 
  
  return(list(gAUCs    = gNewAUCs,
              gAUCi    = gIndividualAUCs,
              plotData = rocAll,
              AUCdata  = retAUC))
}

###### Start
iData <- read.table(paste0(outdir,
                            '/Rtables/k4me3Table.forR.tab'),
                    header=TRUE)

## Pick subset of hist. mods. 
iData <- iData[!(iData$TSS & iData$HS) & iData$cs %in% paste0('chr',1:19),]
iData$qOK <- TRUE

## remove duplicates
## Pick subset of hist. mods. 

hmList <- c("H3K4me3" ,"H3K9ac"  ,"H3K36me3","H3K4me2","H3K4me1",
            "H3K27ac" ,"H4ac5"   ,"H4K20me3","H4K12ac","H4K8ac",
            "H3K79me1","H3K27me3","H3K79me3","H3K4ac" ,"H3K27me1",
            "H3K9me3" ,"H3K9me2" )

## Set types
iData$name                         <- 'Other'
iData$name[iData$TSS  & !iData$HS] <- 'TSS'
iData$name[!iData$TSS & iData$HS]  <- 'HS'

col2Use <- names(iData) %in% c(hmList,'name')
row2Use <- iData$cs %in% paste0("chr",1:19) & ((iData$HS & !iData$TSS) | (iData$TSS & !iData$HS))

iData <- iData[row2Use,col2Use] 
iData <- iData %>% select(.,c(hmList,'name'))
iData$noZeros <- apply(iData[,1:length(hmList)],1,FUN=function(x){sum(x<=0)}) < 1
iOK <- iData[iData$noZeros,1:length(hmList)]


## Omit big outliers
iData$qOK <- TRUE
for (i in hmList){
  iData$qOK <- iData$qOK & (iData[[i]] < quantile(iData[[i]],.99))
}

hmData <- iData[iData$qOK,]

## To TPM
for (i in hmList){
  hmData[[i]] <-  toTPM(hmData[[i]])
}

hmPCAList <- c("H3K9ac"  ,"H3K4me2","H3K4me1",
               "H3K27ac" ,"H4ac5"   ,"H4K20me3","H4K12ac","H4K8ac",
               "H3K79me1","H3K27me3","H3K79me3","H3K4ac" ,"H3K27me1",
               "H3K9me3" ,"H3K9me2" )

rData <- hmData %>% select(.,c(hmPCAList,'name'))
rDim  <- dim(rData)[2]-1

set.seed(18) ### Means "must get rich" in Chinese
objPCA          <- prcomp(rData[,1:rDim],scale=TRUE,center = TRUE)
pcByType        <- cbind(as.data.frame(objPCA$x),rData$name)
names(pcByType) <- c(names(pcByType)[1:rDim],'type')

pcaIndices <- get_pca_ind(objPCA)

pcVariance      <- apply(objPCA$x, 2, var)  
pcVarPercent    <- pcVariance / sum(pcVariance)
pcVarCumulative <- data.frame(var=cumsum(pcVarPercent),PC=1:rDim)

## get PCA dets
pcVarz           <- get_pca_var(objPCA)

dfContrib        <- as.data.frame(pcVarz$contrib)
dfContrib$type   <- rownames(dfContrib)
names(dfContrib) <- c(paste0('PC',1:rDim),'type')

#########
dfPCA           <- cbind(as.data.frame(pcaIndices$coord),name=rData$name)

######### DRAW
noMarg          <- theme(plot.margin=unit(c(0,0,0,0),'cm'))

gROCs <- makeROCsPlot(histmods = hmList,
                      dataHM = hmData, 
                      dataPCA = dfPCA)

gROCs$AUCdata$detail <- factor(gROCs$AUCdata$detail,levels=rev(gROCs$AUCdata$detail))
gROCs$AUCdata$type   <- factor(gROCs$AUCdata$type,levels=rev(levels(gROCs$AUCdata$type)))

gAUC <- ggplot(gROCs$AUCdata, aes(x=AUC, y=detail, label=round(AUC,2))) + 
  geom_vline(xintercept=0.5,lwd=.2,lty='dashed') +
  geom_segment(data=gROCs$AUCdata,aes(xend = 0.5,
                                      yend = detail),
               color = "black", lwd=.2)+
  geom_point(stat='identity',
             shape=21, size=3,
             lwd=.2,
             aes(fill=AUC))+
  geom_vline(xintercept=c(0.25,0.75),lty='dotted',lwd=.1)+
  geom_text(color="black", size=7*5/14, aes(x=AUC+0.08)) +
  scale_fill_gradient2(low='white',high='blue',midpoint=0.55) +
  theme(legend.position='none',
        strip.background = element_blank(),
        strip.text=element_blank()) +
  scale_x_continuous(breaks=seq(0.5,1,0.1)) +
  xlab('Hotspots AUC') +
  ylab('Property') +
  facet_wrap(~type,ncol=2,nrow = 2,scales='free_y') +
  coord_cartesian(xlim=c(0.46,1))

gCumVar <- ggplot(pcVarCumulative,aes(x=PC,y=var*100)) +
  geom_line(lwd=.2) +
  geom_point(shape=21,fill='orange',size=0.9,lwd=.2) +
  scale_y_continuous(breaks=seq(0,100,25)) +
  scale_x_continuous(breaks=seq(1,17,3)) +
  xlab('Principal Component') +
  ylab('Cumulative variance explained (%)') +
  coord_cartesian(xlim=c(1,rDim))

gS11BC <- ggarrange(gAUC,gCumVar,
                    ncol=2,nrow=1,
                    widths=c(2,1),
                    align='h',
                    labels     = c('b','c'),
                    vjust      = 1,
                    hjust      = 0,
                    font.label = list(size = 8,face='bold'))

gS11 <- ggarrange(gROCs$gAUCi,gS11BC,
                  ncol=1,nrow=2,
                  labels     = c('a',''),
                  heights=c(6,4),
                  vjust      = 1,
                  hjust      = 0,
                  font.label = list(size = 8,face='bold'))

### Plot Supp Fig 11
png(file = getIMGname(saveLocation = imgdir,
                      fname = 'Figure_S11',
                      type = 'PNG'),
    width=6,height=7,
    units='in',res=300)
print(gS11)
dev.off()

pdf(file = getIMGname(saveLocation = imgdir,
                      fname = 'Figure_S11',
                      type = 'PDF'),
    width=6,height=7)
print(gS11)
dev.off()

## For figure 4
varData       <- reshape2:::melt.data.frame(dfContrib,id.vars = 'type')
varData$value <- round(varData$value)
varData$value[varData$value<100/rDim] <- NA

varData$lbl   <- as.character(varData$value)
varData$lbl[is.na(varData$lbl)] <- 'ns'

varData$pc    <- as.numeric(substr(varData$variable,3,4))

varData$type <- factor(varData$type,levels=hmList)

gContribs <- ggplot(varData,aes(y=type,x=pc,fill=value)) +
  geom_tile(na.rm = TRUE) +
  geom_vline(xintercept=seq(1.5,rDim+0.5,1),
             lwd=.05,
             lty='dashed')+
  scale_fill_gradient2(midpoint=1,
                       high='royalblue2',
                       low='white',
                       na.value = 'white')+
  geom_text(size=7*5/14,
            aes(color=is.na(value),
                label=lbl)) +
  scale_x_continuous(breaks=seq(1,rDim,2)) +
  theme(legend.position='none') +
  scale_color_manual(values=c('black','grey90')) +
  coord_cartesian(xlim=c(1-0.5,rDim+0.5),expand = FALSE) +
  xlab('Principal component') +
  ylab('Contribution to PC (%)')



