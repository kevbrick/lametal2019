### Camerini Group functions
### KB June 29 2018
#############################

scriptFolder    <- paste0(Sys.getenv('SHARE'),'/code/R/')
imgOutputFolder <- paste0(Sys.getenv('RFIGS'),'/')

## Fix double slashes in foldernames
scriptFolder    <- gsub('//','/',scriptFolder)
imgOutputFolder <- gsub('//','/',imgOutputFolder)

## Load bioconductor (## NOT USED CURRENTLY)
source('http://bioconductor.org/biocLite.R')

## Load libraries
library(plyr)
library(dplyr)
library(data.table)
library(extrafont)
library(factoextra)
library(ggplot2)
library(ggpmisc)
library(ggpubr)
library(grid)
library(gridExtra)
library(numform) ## Format numbers nicely
library(preprocessCore)
library(reshape2)
library(tictoc)

############################# 
### standardizeMNSD
### KB June 29 2018
### Normlize the rows of a data frame by mean & sd
### ARGS: 
# z     numeric data frame to normalize
### OUTPUTS:
# Normalized data frame
standardizeMNSD <- function(z) {
  if (is.null(dim(z))){
    rv <- (z-mean(z))/sd(z)
  }else{
    rowmean <- apply(z, 1, mean)
    rowsd <- apply(z, 1, sd)  
    
    rv <- sweep(z, 1, rowmean,"-")  #subtracting mean expression
    rv <- sweep(rv, 1, rowsd, "/")  # dividing by standard deviation
  }
  return(rv)
}

############################# 
### standardizeMedMad
### KB June 29 2018
### Normlize the rows of a data frame by median & mean absolute deviation
### ARGS: 
# z     numeric data frame to normalize
### OUTPUTS:
# Normalized data frame
standardizeMedMad <- function(z) {
  rowmed <- apply(z, 1, median)
  rowmad <- apply(z, 1, mad)  # median absolute deviation
  
  rv <- sweep(z, 1, rowmed,"-")  #subtracting median expression
  rv <- sweep(rv, 1, rowmad, "/")  # dividing by median absolute deviation
  return(rv)
}

############################# 
### standardize0to1
### KB June 29 2018
### Normlize the rows of a data frame from 0 to 1
### ARGS: 
# z     numeric data frame to normalize
### OUTPUTS:
# Normalized data frame
standardize0to1 <- function(z) {
  ## For vector
  if (is.null(dim(z))){
    rv <- (z-min(z))/(max(z)-min(z))
  }else{ ## For DF / matrix
    rowmin <- apply(z, 1, min)
    rowmax <- apply(z, 1, max)  
    
    rv <- sweep(z, 1, rowmin,"-")  
    rv <- sweep(rv, 1, rowmax-rowmin, "/")  
  }
  return(rv)
}

############################# 
### theme7point
### KB June 29 2018
### Set the default theme to 7-point font
### suitable for publication
### ARGS: none
### OUTPUTS: Sets the theme
theme7point <- function() {
  theme_base_size <- 7
  theme_set(theme_bw(base_size = theme_base_size) %+replace% 
              theme(axis.text = element_text(size=theme_base_size,
                                             color='black'),
                    panel.grid = element_blank(),
                    panel.border=element_blank(),
                    axis.line=element_line(size=.2),
                    axis.ticks = element_line(size=.2)))
}

############################# 
### theme7point
### KB June 29 2018
### Set the default theme to 7-point font
### suitable for publication
### ARGS: none
### OUTPUTS: Sets the theme
theme28point <- function() {
  theme_base_size <- 28
  theme_set(theme_bw(base_size = theme_base_size) %+replace% 
              theme(axis.text = element_text(size=theme_base_size,
                                             color='black'),
                    panel.grid = element_blank(),
                    panel.border=element_blank(),
                    axis.line=element_line(size=.8),
                    axis.ticks = element_line(size=.8)))
}
############################# 
### getImgName
### KB June 29 2018
### Create a systematic image name with date and user
### ARGS: 
# fname         Output file name stem
# type          PNG/PDF/GIF/TIF (default: PDF)
# saveLocation  Folder (default uses system RFIGS path)
### OUTPUTS: File name
getIMGname <- function(fname=NULL,
                       type="PDF", ## 
                       saveLocation=NULL) {
  
  systemUserName <- as.character(Sys.info()[7])
  
  if (is.null(saveLocation)){
    saveLocation <- imgOutputFolder
  }
  
  if (type == 'PDF'){return(paste0(saveLocation,format(Sys.time(), "%y%m%d"),"_",systemUserName,"_",fname,".pdf"))}
  if (type == 'PNG'){return(paste0(saveLocation,format(Sys.time(), "%y%m%d"),"_",systemUserName,"_",fname,".png"))}
  if (type == 'GIF'){return(paste0(saveLocation,format(Sys.time(), "%y%m%d"),"_",systemUserName,"_",fname,".gif"))}
  if (type == 'TIF'){return(paste0(saveLocation,format(Sys.time(), "%y%m%d"),"_",systemUserName,"_",fname,".tiff"))}
}

############################# 
### ggCorMat
### KB July 03 2018
### Generate a ggplot correlation matrix
### ARGS: 
# mCC           correlation matrix
# newOrd1       Order of field 1
# newOrd2       Order of field 2
# flipIt        Flip matrix
# numOFF        Do NOT print CCs
# tileFontScale Scale tile text font (defauilt=1)
# varColor      Use different colors for CCs > 0.5 
# tileFontSize Like it says ... default = 7
### RETURNS: ggplot Geom
ggCorMat <- function(mCC,
                     newOrd1=NULL,
                     newOrd2=NULL,
                     flipIt=FALSE,
                     numOFF=FALSE,
                     scalePC=1,
                     tileFontScale=1,
                     varColor=FALSE,
                     tileFontSize=7,
                     keepLeadingZeros=FALSE,
                     decimalPlaces=1,
                     asPercentage=FALSE,
                     noDiagonal=FALSE,
                     xTilt=45,
                     yOnRight=FALSE){
  library(ggplot2)
  library(reshape2)
  
  ## Set default order
  if (is.null(newOrd1)){
    if (is.data.frame(mCC)){
      newOrd1 <- names(mCC)
    }
    
    if (is.matrix(mCC)){
      newOrd1 <- colnames(mCC)
    }
  }
  
  if (is.null(newOrd2)){
    if (is.data.frame(mCC)){
      newOrd2 <- names(mCC)
    }
    
    if (is.matrix(mCC)){
      newOrd2 <- rownames(mCC)
    }
  }
  
  mCC <- round(mCC,2)
  melted_cormat <- melt(mCC)
  
  # Get lower triangle of the correlation matrix
  get_upper_tri<-function(cormat){
    cormat[upper.tri(cormat)] <- NA
    return(cormat)
  }
  
  # Get upper triangle of the correlation matrix
  get_lower_tri <- function(cormat){
    cormat[lower.tri(cormat)]<- NA
    return(cormat)
  }
  
  if (flipIt){
    lower_tri     <- get_lower_tri(mCC)
    melted_cormat <- melt(lower_tri)
  }else{
    lower_tri     <- get_upper_tri(mCC)
    melted_cormat <- melt(lower_tri)
  }
  
  
  if (length(newOrd1)>0){
    melted_cormat$Var1 <- factor(melted_cormat$Var1,newOrd1)
    if (length(newOrd2)>0){
      melted_cormat$Var2 <- factor(melted_cormat$Var2,newOrd2)
    }else{
      melted_cormat$Var2 <- factor(melted_cormat$Var2,newOrd1)
    }
  }
  
  ## Pad with required number of zeros
  if (asPercentage){
    melted_cormat$numz <- f_num(melted_cormat$value*100, digits=decimalPlaces)      
  }else{
    melted_cormat$numz <- f_num(melted_cormat$value ,digits=decimalPlaces)  
  }
  
  ## Remove leading zeros
  if (keepLeadingZeros){
    melted_cormat$numz <- f_pad_zero(melted_cormat$numz,
                                     width=(decimalPlaces+2))      
  }
  
  colScale <- colorRampPalette(c("green", "blue", "black", "#FF8811", "red"))
  
  if (varColor){
    melted_cormat$txtColor <- 0
    melted_cormat$txtColor[melted_cormat$value > 0.5] <- 1
  }else{
    melted_cormat$txtColor <- 0
  }
  
  if (asPercentage){
    myLbl <- 'Overlap (%)'
  }else{
    myLbl <- substitute(paste('Spearman ', R^"2"))
  }
  
  if (noDiagonal){
    melted_cormat <- melted_cormat[melted_cormat$Var1 != melted_cormat$Var2 &
                                     melted_cormat$Var2 != newOrd1[1] & 
                                     melted_cormat$Var1 != newOrd2[length(newOrd2)],]
  }
  
  if (xTilt == 45){
    xVjust <- 1 
    xHjust <- 1 
  }
  
  if (xTilt == 0){
    xVjust <- .5 
    xHjust <- .5 
  }
  
  ggheatmap <-  ggplot(melted_cormat, aes(Var2, Var1, fill = value)) +
    geom_tile(color = "white",size=.3) +
    scale_fill_gradient2(low = "grey90", mid = 'orange', high = "orangered3", 
                         midpoint = 0.5, limit = c(0,1), space = "Lab", 
                         name=myLbl,na.value='white') +
    theme(axis.text.x = element_text(vjust = xVjust, 
                                     hjust = xHjust,
                                     angle = xTilt,
                                     size = tileFontSize),
          axis.text.y = element_text(size = tileFontSize)) 
  
  coord_fixed()
  #    theme(axis.text = element_text(size = tileFontSize)) +
  
  gCInit <- ggheatmap + 
    theme(
      axis.title.x         = element_blank(),
      axis.title.y         = element_blank(),
      panel.grid.major     = element_blank(),
      panel.border         = element_blank(),
      panel.background     = element_blank(),      
      axis.ticks           = element_blank(),
      legend.position      = 'none',
      legend.justification = c(1,1),
      legend.text          = element_text(size = tileFontSize),
      legend.title         = element_text(size = tileFontSize),
      legend.key.size      = unit(.5,'cm'),
      legend.direction     = "vertical")
  #legend.key.size = unit(1*scalePC,'cm'),
  #legend.key.width = unit(0.6*scalePC,'cm'))
  
  if (numOFF){
    gCC <- gCInit
  }else{
    ## NOTE: Font sizes in geom_text are 14/5 X bigger than in theme 
    ## NO idea why !!!
    
    gCC <- gCInit + 
      geom_text(aes(Var2, Var1, label = numz, color = txtColor), 
                size = 5/14*tileFontSize,
                show_guide = FALSE) + 
      scale_color_gradient2(low = "black", mid = 'black', high = "white", 
                            midpoint = 0.1, limit = c(0,1))
  }
  
  if (yOnRight){
    gCC <- gCC + scale_y_discrete(position='right')
  }
  
  return(gCC)
}

############################# 
### toTPM
### KB July 03 2018
### Convert a vector of strengths into a Tags (or Fragments) per million value
### ARGS: 
# x           vector of values
# noNeg       remove negative values by adding the minimum FPM to all 
## RETURNS: Normalized vector in T(F)PM
toTPM <- function (x,noNeg=FALSE){
  #x <- x+abs(min(x))+1;
  v <- (x/sum(x)*1000000)
  if (noNeg & min(v) < 0){
    v <- v - min(v) + 1
  }
  return(v)
}

############################# 
### convertToQuantiles
### KB July 03 2018
### Convert a vector of values into evenly size quantiles
### ARGS: 
# pData       vector of values
# nQ          number of quantiles (default = 5)
# plotMe      return a diagnostic plot
# revLbl      reverse order (default == smallest = 1)
# labelMaxMin change label of max and min bins to "max" and "min"
# maxLabel    max Label
# minLabel    min Label
## RETURNS: a factorized vector of quantiles

convertToQuantiles <- function (pData,nQ=5,plotMe=FALSE,revLbl=FALSE,
                                labelMaxMin=FALSE,
                                maxLabel=NULL,
                                minLabel=NULL){
  library(lsr)
  library(ggplot2)
  
  qvec <- quantile(jitter(pData), 
                   probs=seq(0,1,length=(nQ+1)), 
                   type=1,
                   na.rm=T)
  
  # Old alternate (KB)
  #asdc <- function(q){max(1,which(qvec<=q))}
  #qvals = apply(matrix(pData,nrow=1),2,asdc)
  
  if (revLbl){
    qvals <- cut(pData, breaks = qvec,include.lowest = TRUE, labels = nQ:1)
    qvals <- quantileCut(jitter(pData,0.001),nQ,labels=nQ:1)
  }else{
    qvals <- cut(pData, breaks = qvec,include.lowest = TRUE, labels = 1:nQ)
    qvals <- quantileCut(jitter(pData,0.001),nQ,labels=1:nQ)
  }
  
  
  
  if (plotMe){
    df <- data.frame("pD"=pData,"Q"=qvals)
    g <- ggplot(df,aes(x=Q,y=..count..,color=Q)) + geom_bar() 
    + theme(text=element_text(size=22)) 
    + xlab("Quantile") 
    + ylab("Count")
    print(g)
  }
  
  if (labelMaxMin){
    if (is.null(maxLabel)){maxLabel = 'Max'}
    if (is.null(minLabel)){minLabel = 'Min'}
    
    qV            <- qvals
    qV            <- factor(qV,levels=c(minLabel,2:(nQ-1),maxLabel))
    
    qV[qvals==nQ] <- maxLabel
    qV[qvals==1]  <- minLabel
    qvals <- qV
  }
  
  return(qvals)
}

############################# 
### getLegend
### KB July 10 2018
### get the legend of a ggplot
### ARGS: 
# g   grob
## RETURNS: a legend as a grob
getGGLegend<-function(g){
  tmp <- ggplot_gtable(ggplot_build(g))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

############################# 
### saveImage
### KB July 10 2018
### Save a plot as a PNG and PDF
### ARGS: 
# imgName      name of image
# img2Output   plot
# nH           height
# nW           width
## RETURNS: Nothing
saveImage <- function(imgName, img2Output, nH, nW){
  outPNG <- getIMGname(imgName,'PNG')
  outPDF <- getIMGname(imgName,'PDF')
  
  cat ('************* ERROR ***************\n')
  cat ("saveImage function NOT WORKING !!!!\n")
  cat ("***********************************\n")
  #### Currently NOT working
  return()
  
  ## Closes active graphics devices
  graphics.off()
  
  png(filename = outPNG, width = nH, height = nW, units='in', res=300)
  img2Output
  dev.off()
  
  pdf(file = outPDF, width = nH, height = nW)
  img2Output
  dev.off()
  
  print(paste0('Saved as ',outPNG))
  print(paste0('Saved as ',outPDF))
}

############################# 
### fancy_scientific
### KB July 10 2018
### Change number to scientific notation format:
### Most useful for scale_x[y]_log10(labels=fancy_scientific)
### ARGS: 
# l   number
## RETURNS: Formatted number as expression
fancy_scientific <- function(l) {
  # turn in to character string in scientific notation
  l <- format(l, scientific = TRUE)
  # quote the part before the exponent to keep all the digits
  l <- gsub("^(.*)e", "e", l)
  # remove +s
  l <- gsub("\\+", "", l)
  # turn the 'e+' into plotmath format
  l <- gsub("e", "10^", l)
  # return this as an expression
  parse(text=l)
}

############################# 
### addThemeToStatCompare
### KB July 10 2018
### Fix font size and line width in stat_compare
### ARGS: 
# f   ggplot figure
## RETURNS: ggplot figure
addThemeToStatCompare <- function(f,fLineSize=0.2,fFontSize=7) {
  f$layers[[which_layers(f, "GeomSignif")]]$aes_params$size <- fLineSize
  f$layers[[which_layers(f, "GeomSignif")]]$aes_params$textsize <- fFontSize*5/14
  return(f)
}

############################# 
### darkenColor
### KB July 17 2018
### Darken a color
### ARGS: 
# color   color
# X       how much darker? (default=1.4x)
## RETURNS: color
darkenColor <- function(color, X=1.4){
  col <- col2rgb(color)
  col <- col/X
  col <- rgb(t(col), maxColorValue=255)
  col
}

############################# 
### lightenColor
### KB July 17 2018
### Lighten a color
### ARGS: 
# color   color
# X       how much lighter? (default=1.4x)
## RETURNS: color
lightenColor <- function(color, X=1.4){
  col <- col2rgb(color)
  col <- col*X
  col <- rgb(t(col), maxColorValue=255)
  col
}

############################# 
### plotHSstrength
### KB July 23 2018
### Make a density scatterplot of hotspot strength in two samples
### ARGS: 
# sA      vector of values for x-axis
# sB      vector of values for y-axis
# nameA   Name for x-axis
# nameB   Name for y-axis
# 
## RETURNS: geom
plotHSstrength <- function (sA,sB,nameA,nameB,
                            noGradient = FALSE,
                            sOK        = NULL,
                            lblAsIs    = FALSE,
                            topLeft    = FALSE,
                            xLim       = 0,
                            yLim       = 0,
                            xLimits    = NULL,
                            yLimits    = NULL,
                            stdStyle   = TRUE,
                            pointSz    = 2,
                            txtSz      = 7,
                            noLeg      = FALSE){
  
  sDataInit <- data.frame('tpmA'=toTPM(sA),'tpmB'=toTPM(sB))  
  
  if (length(sOK) > 0){
    sData <- sDataInit[sOK>0,]
  }else{
    sData <- sDataInit
  }
  
  ## Get max / min coords
  xMax <- max(sData$tpmA)
  xMin <- min(sData$tpmA)
  
  yMax <- max(sData$tpmB)
  yMin <- min(sData$tpmB)
  
  ## Set limits if we're using them
  if (xLim){
    xMax <- xLim
  }
  
  if (yLim){
    yMax <- yLim
  }
  
  if (!is.null(xLimits)){
    xMin <- xLimits[1]
    xMax <- xLimits[2]
  }
  
  if (!is.null(yLimits)){
    yMin <- yLimits[1]
    yMax <- yLimits[2]
  }
  
  sCC <- round(cor(sData[,1:2],method='spearman')[1,2]^2,2)
  
  if (lblAsIs){
    xLblName <- nameA
    yLblName <- nameB
  }else{
    xLblName <- paste0(nameA,' ','(FPM)')
    yLblName <- paste0(nameB,' ','(FPM)')
  }
  
  okData <- sData[sData$tpmA>xMin & sData$tpmA<xMax & 
                    sData$tpmB>yMin & sData$tpmB<yMax,]
  
  if (stdStyle){
    gInit <- ggplot(data=okData,aes(x=tpmA,y=tpmB)) + 
      geom_point(color='grey50',size=pointSz) +
      scale_x_log10(labels=fancy_scientific) + 
      scale_y_log10(labels=fancy_scientific) + 
      geom_smooth(method=lm,linetype=2,colour="NA",se=F) + 
      guides(alpha="none") + 
      annotation_logticks(sides='lb',
                          size=.2,
                          short=unit(0.050,'cm'),
                          mid=unit(0.075,'cm'),
                          long=unit(0.100,'cm')) +
      coord_cartesian(xlim=c(xMin,max(xMax,yMax)),
                      ylim=c(xMin,max(xMax,yMax))) +
      xlab(xLblName) + 
      ylab(yLblName) + 
      theme(legend.position=c(1,0),
            legend.justification=c(1,0)) 
  }else{
    gInit <- ggplot(data=okData,aes(x=tpmA,y=tpmB)) + 
      geom_point(color='grey50',size=pointSz) +
      scale_x_log10(labels=fancy_scientific) + 
      scale_y_log10(labels=fancy_scientific) + 
      geom_smooth(method=lm,linetype=2,colour="NA",se=F) + 
      guides(alpha="none") + 
      annotation_logticks(sides='lb') +
      coord_cartesian(xlim=c(xMin,max(xMax,yMax)),
                      ylim=c(xMin,max(xMax,yMax))) +
      theme_MF() + 
      xlab(xLblName) + 
      ylab(yLblName) + 
      theme(plot.title=element_text(size=20),
            legend.position=c(1,0),
            legend.justification=c(1,0),
            legend.title=element_text(size=15)) 
  }
  #ggtitle(substitute(paste('Spearman ', R^"2"," = ", cc ,sep=''),list(cc = sCC)))
  
  if (noGradient){
    g <- gInit
  }else{
    g <- gInit + stat_density2d(aes(fill=..level..,
                                    alpha=..level..),
                                geom='polygon',
                                colour='NA')  
    #scale_fill_continuous(low="grey30",high="red") 
  }  
  
  g$labels$fill <- "Hotspot density"
  
  myLbl <- substitute(paste('Spearman ', R^"2"," = ", cc ,sep=''),list(cc = sCC));
  myLbl <- paste("R^2 == ", sCC)
  labelSize <- 10
  if (topLeft){
    if (stdStyle){
      
      thisTheme <- theme_get()
      txtSz <- thisTheme$axis.text$size
      
      gg <- g + annotate("text", 
                         label = myLbl, 
                         x = -Inf, y = Inf, 
                         colour = "black", hjust=0,
                         parse=TRUE,
                         size = 5/14*txtSz) +
        theme(axis.line = element_line(color='black'),
              legend.background=element_blank())
    }else{
      gg <- g + annotate("text", 
                         label = myLbl, 
                         x = -Inf, y = Inf, 
                         size = labelSize, 
                         colour = "black", hjust=0,
                         parse=TRUE) +
        theme(axis.line = element_line(color='black'),
              legend.background=element_blank())
    }
  }else{
    if (stdStyle){
      gg <- g + annotate("text", 
                         label = myLbl, 
                         x = yMin, y = yMax-(yMax*.1),  
                         colour = "black", hjust=0,
                         parse=TRUE) +
        theme(axis.line = element_line(size=.5),
              legend.background=element_blank())
    }else{
      gg <- g + annotate("text", 
                         label = myLbl, 
                         x = yMin, y = yMax-(yMax*.1),  
                         size = labelSize, 
                         colour = "black", hjust=0,
                         parse=TRUE) +
        theme(axis.line = element_line(size=.5),
              legend.background=element_blank())
    }
  }
  
  if (noLeg){
    gg <- gg + theme(legend.position='none')
  }
  
  return(gg)
}

############################# 
### makeROC
### KB Oct 13 2018
### Make an ROC curve
### ARGS: 
# score   ranked scores
# type    corresponding TRUE/FALSE positives
# 
## RETURNS: geom, AUC
makeROC <- function(score,type,tryFR=FALSE,title=NULL){
  
  library(pROC)
  
  getROC <- function(sc,t){
    dfInit <- data.frame(score=sc,
                         type=t)
    
    dfS     <- dfInit[order(sc, decreasing = TRUE),]
    dfS$TPR <- cumsum(dfS$type) / sum(dfS$type)
    dfS$FPR <- cumsum(!dfS$type) / sum(!dfS$type)
    
    roc_obj <- roc(dfS$type,dfS$score)
    AUC     <- auc(roc_obj)
    manualAUC <- sum(dfS$TPR)/length(dfS$TPR)
    
    return(list(df = dfS,
                ro = roc_obj,
                AUC1 = AUC,
                AUCm = manualAUC))
  }
  
  if (tryFR){
    fROC <- getROC(score,type)
    rROC <- getROC(score*-1,type)
    
    if (fROC$AUCm >= rROC$AUCm){
      myROC <- fROC
    }else{
      myROC <- rROC
    }
    
  }else{
    myROC <- getROC(score,type)
  }
  
  dfS <- myROC$df
  
  theme7point()
  gROC <- ggplot(dfS,aes(x=FPR,y=TPR)) + 
    geom_line(data=data.frame(xV=c(0,1),
                              yV=c(0,1)),
              aes(x=xV,y=yV),
              color='grey30') + 
    geom_line(color='red') + 
    geom_text(x=0.6,y=0.4,
              label=paste0(title,"AUC = ",round(myROC$AUCm,2)),
              size=7*5/14,
              check_overlap = TRUE)
  
  return(list(data=dfS,
              ROC=gROC,
              AUC=myROC$AUC1,
              manualAUC=myROC$AUCm))
}

####################################################################################### 
### ndrDistrib
### ADDED 01/02/2019: KB
### Generate a distribution of ndr locations (initially for modelling Spo11 cutting)
### Sequential penalties are incurred for each nucleosome traversed
### ARGS: 
# nNuc    number of nucleosome to model either side of center
# nSD     standard deviation (affects shape of distrib)
# pNext   Probability decrease for each adjacent nucleosome
#         (i.e. 0.1 -> the +1 NDR is 10% less likely used than the 0 NDR)
#         (            the +2 NDR is 10% less likely used than the +1 NDR)
#         (            the +3 NDR is 10% less likely used than the +2 NDR)
#         (            .... etc ... )
# 
## RETURNS: dataframe of values (for large N, can get very big)
ndrDistrib <- function(nNuc=10,nSD=20,pNext=.1){
  nZero   <- 100000
  dCenter <- rnorm(nZero,0,sd = nSD)
  
  posNuc  <- 0
  probNuc <- 1
  
  for (n in 2:nNuc){
    posNuc  <- posNuc + 150
    probNuc <- probNuc * pNext
    
    dCenter <- c(dCenter,
                 rnorm(nZero*probNuc, posNuc   , sd = nSD),
                 rnorm(nZero*probNuc, -1*posNuc, sd = nSD))
  }
  
  return(data.frame(N=dCenter))
}

####################################################################################### 
### nucDistrib
### ADDED 01/02/2019: KB
### Generate a distribution of nucleosome locations 
### Sequential penalties are incurred for each nucleosome traversed
### NOTE: this was designed to model the positions excluded by ndrDistrib
### ARGS: 
# nNuc    number of nucleosome to model either side of central nucleosome
# nSD     standard deviation (affects shape of distrib)
# pNext   Probability decrease for each adjacent nucleosome
#         (i.e. 0.1 -> +1 is 10% less likely used than 0 )
#         (            +2 is 10% less likely used than +1 )
#         (            +3 is 10% less likely used than +2 )
#         (            .... etc ... )
# 
## RETURNS: dataframe of values (for large N, can get very big)
nucDistrib <- function(nNuc=10,nSD=20,pNext=.1){
  nZero   <- 100000
  dCenter <- c(rnorm(nZero,-75,sd = nSD),
               rnorm(nZero,75,sd = nSD))
  
  posNuc  <- 75
  probNuc <- 1
  
  for (n in 2:nNuc){
    posNuc  <- posNuc + 150
    probNuc <- probNuc * pNext
    
    dCenter <- c(dCenter,
                 rnorm(nZero*probNuc, posNuc   , sd = nSD),
                 rnorm(nZero*probNuc, -1*posNuc, sd = nSD))
  }
  
  return(data.frame(N=dCenter))
}

####################################################################################### 
### getDensityDF
### ADDED 01/02/2019: KB
### Generate data for a density-smoothed histogram
### Faster than running geom_density every time (esp. for large lists)
### ARGS: 
# dd      data frame of values (required cols = from, to, strand)
# nFr     start density from
# nTo     end density at
#
## RETURNS: dataframe of values
getDensityDF <- function(dd,nFr,nTo){
  dd$mid <- round((dd$from + dd$to) / 2)
  zDN <- density(dd$mid[dd$strand == -1],from=nFr,to=nTo)
  zDP <- density(dd$mid[dd$strand ==  1],from=nFr,to=nTo)
  
  zBoth <- data.frame(x=c(zDP$x,zDN$x),
                      y=c(zDP$y,zDN$y),
                      strand=c(rep(1,length(zDP$x)), 
                               rep(-1,length(zDN$x))))
  return(zBoth)
}

####################################################################################### 
### getRMSE
### ADDED 01/02/2019: KB
### Calculate the Root Mean Squared Error between two vectors
### ARGS: 
# m      first vector
# o      second vector
## RETURNS: RMSE value
getRMSE = function(m, o){
  NAok <- !is.na(m*o)
  NAok <- NAok & (m*o > 0)
  retRMSE <- sqrt(mean((m[NAok] - o[NAok])^2))
  return(retRMSE)
}

####################################################################################### 
### getChromSize
### ADDED 01/24/2019: KB
### Get chromosome size
### ARGS: 
# myGenome  genome
# myCS      chromosome
## Get chromosome sizes
getChromSize <- function(myGenome,
                         myCS,
                         noSexCS=FALSE,
                         noX=FALSE,
                         noY=FALSE,
                         noM=TRUE){
  
  '%!in%' <- function(x,y)!('%in%'(x,y))
  
  allCS <- c(paste0('chr',c(1:50)),'chrX','chrY','chrM')
  
  if (noSexCS){allCS <- allCS[allCS %!in% c('chrX','chrY','chrZ')]}
  if (noM){allCS <- allCS[allCS %!in% 'chrM']}
  if (noX){allCS <- allCS[allCS %!in% 'chrX']}
  if (noY){allCS <- allCS[allCS %!in% 'chrY']}
  
  genomeSz         <- read.table(paste0(Sys.getenv('GENOMES'),'/',myGenome,'/genome.fa.fai'),header=FALSE)
  names(genomeSz)  <- c('cs','size','cum','xa','xb')
  
  genomeSz$size    <- as.numeric(genomeSz$size)
  
  genomeSz$cs      <- factor(genomeSz$cs,levels=allCS)
  
  genomeSz         <- genomeSz[!is.na(genomeSz$cs),]
  
  genomeSz         <- genomeSz[order(genomeSz$cs),]
  
  genomeSz$min     <- c(0,cumsum(genomeSz$size[1:(length(genomeSz$size)-1)]))+1
  genomeSz$max     <- cumsum(genomeSz$size)
  
  genomeSz <- genomeSz[,c('cs','size','min','max')]
  
  ## Remove unused chromosomes
  z <- plyr:::join(data.frame(cs=levels(genomeSz$cs)),genomeSz,by='cs')  
  genomeSz$cs <- factor(genomeSz$cs, levels=z$cs[!is.na(z$size)])
  
  if (myCS %in% c('All','all','ALL')){
    return(genomeSz)
  }else{
    return(genomeSz$size[genomeSz$cs == myCS])
  }
}

## Gives count, mean, standard deviation, standard error of the mean, and confidence interval (default 95%).
##   data: a data frame.
##   measurevar: the name of a column that contains the variable to be summariezed
##   groupvars: a vector containing names of columns that contain grouping variables
##   na.rm: a boolean that indicates whether to ignore NA's
##   conf.interval: the percent range of the confidence interval (default is 95%)
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE, q =.95) {

  library(plyr)
  
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N      = length2  (xx[[col]], na.rm=na.rm),
                     mean   = mean     (xx[[col]], na.rm=na.rm),
                     median = median   (xx[[col]], na.rm=na.rm),
                     pcHi   = as.numeric(quantile (xx[[col]], q,   na.rm=na.rm)),
                     pcLo   = as.numeric(quantile (xx[[col]], 1-q, na.rm=na.rm)),
                     sd     = sd       (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  
  # Rename the "mean" column    
  #datac <- rename(datac, c("mean" = measurevar))
  names(datac)[which(names(datac) == 'mean')] <- measurevar
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}


getStageColors <- function(nAlpha=1,
                           darkLvl=1){
  
  cLE <- alpha(darkenColor('#c87137ff',darkLvl),nAlpha)
  cZY <- alpha(darkenColor('#ffb380ff',darkLvl),nAlpha)
  cEP <- alpha(darkenColor('#ffcc00ff',darkLvl),nAlpha)
  cLP <- alpha(darkenColor('#ffaaeeff',darkLvl),nAlpha)
  cDI <- alpha(darkenColor('#619cffff',darkLvl),nAlpha)
  
  return(c(cLE,cZY,cEP,cLP,cDI))
}


reOrderClusters <- function(xMat,labels2Use){
  ## Renumber clusters to retain pseudo-temporal order
  ## add tiny pseudocount to break ties
    df <- data.frame(k5=xMat$k5,
                   nBest=apply(xMat[,labels2Use],
                               1,
                               function(x){
                                 x <- x +runif(length(x),0.0000001,0.0000004); 
                                 return(which(x==max(x)))
                                 })
                   )
  
  dAgg <- aggregate(df,by = list(df$k5),FUN=mean)
  
  newC  <- 1
  for (dOrd in sort(as.numeric(dAgg$nBest))){
    oldCluster <- dAgg$k5[dAgg$nBest == dOrd]
    dAgg$newCluster[dAgg$nBest == dOrd] <- newC
    xMat$newClusterNumber[xMat$k5==oldCluster] <- newC
    newC = newC + 1
  }
  
  xMat$k5       <- xMat$newClusterNumber
  xMat$cluster  <- xMat$k5
  
  return(xMat)
}


makeCCFigs <- function(pMat,nK,nRand = 10){
  ccs <- data.frame(cc=matrix(data=0,nrow = nK),
                    ccS1=matrix(data=0,nrow = nK),
                    ccS2=matrix(data=0,nrow = nK),
                    k=matrix(data=0,nrow = nK))
  
  nVals <- (nRand*2)+1;
  dCC <- data.frame(r=matrix(data=0,nrow = nK*nVals),
                    k=matrix(data=0,nrow = nK*nVals),
                    type=matrix(data='Unknown',nrow = nK*nVals))
  
  dCC$type <- factor(dCC$type,levels=c('Unknown','R','sGene','sAll'))
  
  pvals <- data.frame(n=1:nK,
                      newClusterNumber=1:nK,
                      p=0,cc=0,pval=1)
  
  for (n in 1:nK){
    cc <- getH3vExpCC(pMat[pMat$k5 == n,],nRand)
    ccs[n,'cc'] <- cc$cc;
    ccs[n,'ccS1'] <- max(cc$geneShuf);
    ccs[n,'ccS2'] <- max(cc$allShuf);
    ccs[n,'k'] <- n;
    d <- data.frame(r=c(cc$cc,cc$geneShuf,cc$allShuf),
                    k=rep(n,nVals),
                    type=c('R',rep('sGene',nRand),rep('sAll',nRand)))
    
    pvals$cc[n] <- cc$cc
#     if (cc$cc > max(cc$geneShuf)){
#       p <- 1/length(cc$geneShuf)
#  #     pvals$p[n] <- paste0('P < ',1/length(cc$geneShuf))
#     }else{
#       p <- sum(cc$geneShuf > cc$cc) / length(cc$geneShuf)
# #      pvals$p[n] <- paste0('P = ',sum(cc$geneShuf > cc$cc) / length(cc$geneShuf))
#     }
    if (cc$cc > 0){
      p <- max(1,sum(cc$geneShuf > cc$cc)) / length(cc$geneShuf)
    }else{
      p <- max(1,sum(cc$geneShuf < cc$cc)) / length(cc$geneShuf)
    }
    
    pvals$pval[n] <- p
    pvals$p[n] <- '(ns)'
    
    if (p <= 0.0001){
      pvals$p[n] <- '***'
    }else{
      if (p <= 0.001){
        pvals$p[n] <- '**'
      }else{    
        if (p <= 0.01){
          pvals$p[n] <- '*'
        }
      }
    }

    dCC[(1+(n-1)*nVals):(nVals*n),] <- d
  }
  
  #geom_point(size=3,shape=21,position=position_jitter(width=.2),aes(fill=factor(k))) + 
  ggRvCluster <- ggplot(dCC,aes(x=type,y=r)) +  
    geom_boxplot(notch=FALSE) + 
    theme(legend.position='none') + 
    xlab('') +
    ylab('Spearman R') 
  
  #geom_point(size=3,shape=21,position=position_jitter(width=.2),aes(fill=factor(type))) + 
  ggRvClusterF <- ggplot(dCC,aes(x=type,y=r)) +  
    geom_boxplot(notch=FALSE) + 
    theme(legend.position='none') + 
    xlab('') +
    ylab('Spearman R')   + facet_grid(.~k) + scale_fill_manual(values=c('orange','grey30','grey70'))
  
  pRet <- list(data   = dCC,
               pvals  = pvals,
               allFig = ggRvCluster,
               cFig   = ggRvClusterF)
}


getH3vExpCC <- function(cMat,nRand = 100){
  
  ## Get CCs 
  trueCC <- cor(c(as.matrix(cMat[,expNames])),
                c(as.matrix(cMat[,stageNames])), 
                method='spearman')
  
  ## By-gene shuffle
  tZero <- tic()
  ccShuf <- matrix(nrow = nRand); 
  for (n in 1:nRand){
    expTSS <- cMat[,expNames]
    stgTSS <- cMat[,stageNames]
    expTSS <- t(apply(expTSS,1,sample))
    stgTSS <- t(apply(stgTSS,1,sample))
    ccShuf[n] <- cor(c(as.matrix(expTSS)),c(as.matrix(stgTSS)), method='spearman')
    
    nRem <- (nRand-n) * ((tic()-tZero)/n)
    #print(paste0('Rep ',n,' ... tRem = ',round(nRem),' s'))
  }
  
  ## Gene-independent shuffle
  nC <- matrix(nrow = nRand); 
  for (n in 1:nRand){
    nC[n] <- cor(c(as.matrix(cMat[,expNames])),sample(c(as.matrix(cMat[,stageNames]))), method='spearman')
  }
  
  lRet <- list(cc=trueCC,
               geneShuf=ccShuf,
               allShuf=nC)
}

drawSlice <- function(sliceLEDZ='',
                      featureSize=NULL,
                      smoothingWindow = 1,
                      winSz=25,
                      RPK=TRUE){
  
  library(png)
  library(data.table)
  
  sliceFolder  <- gsub(x=sliceLEDZ,pattern = 'LE.+tab',replacement='')
  sliceDZ      <- gsub(x=sliceLEDZ,pattern = '^.+LE.',replacement='')
  slicePosA    <- gsub(x=sliceLEDZ,pattern = '^.+.chr',replacement='chr')
  slicePos     <- gsub(x=slicePosA,pattern = '.DZ.tab',replacement='')
  
  grHS         <- read.table(paste0(sliceFolder,'HS.',slicePos,'.bed'))
  names(grHS)  <- c('cs','from','to')
  grHS$mid     <- round((grHS$from+grHS$to)/2)
  
  grHS$minY    <- 1; grHS$maxY <- 2
  grHS$type    <- 'Hotspots'
  
  grTSS        <- read.table(paste0(sliceFolder,'TSS.',slicePos,'.bed'))
  names(grTSS) <- c('cs','from','to')
  grTSS$mid    <- round((grTSS$from+grTSS$to)/2)
  
  grTSS$minY <- 1 ; grTSS$maxY <- 2
  grTSS$type <- 'TSS'
  
  cAnnotation <- rbind(grHS,grTSS)
  
  colfunc<-colorRampPalette(c("gold2","darkorange","firebrick4","grey30","forestgreen","springgreen"))
  colz <- colfunc(6)
  
  #colz<-getStageColors()
  
  aLE <- read.table(paste0(sliceFolder,'LE.',sliceDZ))
  aZY <- read.table(paste0(sliceFolder,'ZY.',sliceDZ))
  aEP <- read.table(paste0(sliceFolder,'EP.',sliceDZ))
  aLP <- read.table(paste0(sliceFolder,'LP.',sliceDZ))
  aDI <- read.table(paste0(sliceFolder,'DI.',sliceDZ))
  
  aLE$stage <- 'LE'
  aZY$stage <- 'ZY'
  aEP$stage <- 'EP'
  aLP$stage <- 'LP'
  aDI$stage <- 'DI'
  
  if (RPK){
    adjFactor <- (smoothingWindow*winSz)/1000
  }else{
    adjFactor <- (smoothingWindow*winSz)
  }
  
  aLE$V3 <- frollmean(aLE$V3,smoothingWindow)/adjFactor
  aZY$V3 <- frollmean(aZY$V3,smoothingWindow)/adjFactor
  aEP$V3 <- frollmean(aEP$V3,smoothingWindow)/adjFactor
  aLP$V3 <- frollmean(aLP$V3,smoothingWindow)/adjFactor
  aDI$V3 <- frollmean(aDI$V3,smoothingWindow)/adjFactor
  
  cAll <- rbind(aLE,aZY,aEP,aLP,aDI)
  
  names(cAll) <- c('cs','pos','cover','stage')
  
  cAll$stage <- factor(cAll$stage,levels=c('LE','ZY','EP','LP','DI'))
  
  gInit <-  
    
    diff <- (max(cAll$pos)-min(cAll$pos))*0.05
  xMin <- (min(cAll$pos)+diff)/1000000
  xMax <- (max(cAll$pos)-diff)/1000000
  
  getMaxVal <- function(x){
    floor(x*10^abs(floor(log10(x))))/10^abs(floor(log10(x)))
  }
  
  ## Zero out BG
  cAll$cover[cAll$cover<=0] <- 0
  cAll$cover[is.na(cAll$cover)] <- 0
  
  ## Add one zero point at either side to ensure that the geom_polygon
  ## has a flat base at zero ... Silly but necessary !! 
  dMinus     <- aggregate(cAll$pos,by=list(cAll$stage),function(x){min(x)-1})
  dPlus      <- aggregate(cAll$pos,by=list(cAll$stage),function(x){max(x)+1})
  addMeMinus <- data.frame(cs=cAll$cs[1],pos=dMinus$x,cover=0,stage=dMinus$Group.1)
  addMePlus  <- data.frame(cs=cAll$cs[1],pos=dPlus$x,cover=0,stage=dPlus$Group.1)
  
  cAll <- rbind(addMeMinus,cAll,addMePlus)
  
  #yMax <- getMaxVal(max(cAll$cover))
  #yMax <- floor(max(cAll$cover)/100)*100
  yMax <- min(pretty(max(cAll$cover)*.9))
  
  gData <- ggplot(cAll,
                  aes(x=pos/1000000,
                      y=cover,
                      fill=stage,
                      color=stage)) + 
    geom_polygon(lwd=0.2) + 
    facet_grid(stage~.) + 
    coord_cartesian(xlim=c(xMin,xMax)) +
    scale_color_manual(values=getStageColors(darkLvl=1.5)) + 
    scale_fill_manual(values=getStageColors()) + 
    xlab('') + 
    ylab('H3K4me3 ChIP-Seq Coverage (A.U.)') + 
    scale_y_continuous(breaks=c(0,yMax),labels=c('',yMax)) +
    theme(legend.position='none',
          strip.background = element_blank(),
          strip.text = element_blank(),
          axis.text.x  = element_blank(),        
          axis.ticks.x = element_blank(),
          axis.line.x = element_blank(),
          plot.margin = unit(c(0,0,-0.3,0),'cm')) + 
    geom_text(aes(x, y, label=lab),
              data=data.frame(x=Inf, y=Inf, 
                              lab=c("Leptonema","Zygonema","Early Pachynema","Late Pachynema","Diplonema"),
                              stage=c("LE","ZY","EP","LP","DI")),
              vjust=1, hjust=1, size=7*5/14)
  
  grHS$n <- 1
  grTSS$n <- 0
  
  grAnnotation <- unique(rbind(grHS,grTSS))
  
  grAnnotation$mid  <- round((grAnnotation$from+grAnnotation$to)/2)
  grAnnotation$from <- grAnnotation$mid-500
  grAnnotation$to   <- grAnnotation$mid+500
  
  gAnnotation <- ggplot(cAll,
                        aes(x=pos/1000000,
                            y=cover)) +
    geom_rect(data=grAnnotation, lwd = 0.2, 
              aes(x=from/1000000,y=to/1000000,
                  xmin=from/1000000,xmax=to/1000000,
                  ymin=1+n,ymax=2+n,
                  fill=type)) +
    coord_cartesian(xlim=c(xMin,xMax),ylim=c(0.95,3.05)) +
    scale_color_manual(values=c(alpha(darkenColor('grey30',1),1),
                                alpha(darkenColor('red',1))),1) + 
    scale_fill_manual(values=c(alpha(darkenColor('grey30',1.3),1),
                               alpha(darkenColor('red',1.3),1))) + 
    xlab(paste0('Position on ',cAll$cs[1],' (Mbp)')) + 
    ylab('') + 
    theme(legend.position='none',
          axis.ticks.y = element_blank(),
          axis.line.y = element_blank(),
          strip.background = element_blank(),
          strip.text = element_blank(),
          plot.margin = unit(c(0,0,0,0),'cm'))+
    scale_y_continuous(breaks=c(1.5,2.5),
                       labels=c('TSS','HS'))
  
  # geom_text(aes(x, y, label=lab),
  #             data=data.frame(x=-Inf, y=1.5, 
  #                             minX=Inf,maxX=Inf,
  #                             minY=Inf,maxY=Inf,
  #                             lab=c("Hotspots","TSS"),
  #                             type=c("Hotspots","TSS")),
  #             vjust=0.5, hjust=1, size=7*5/14)
  
  gCoverage <-ggarrange(gData,gAnnotation, 
                        align = 'v',
                        heights=c(10,2),
                        ncol=1,nrow=2)
  
  return(list(gFig    = gCoverage,
              gTracks = gData,
              gAnno   = gAnnotation)) 
}

nanRows <- function(z) {
  nanRow <- which(apply(dfClust[,stageNames], 1, function(x){sum(is.nan(x))}) > 0)
  return(nanRow)
}

notNanRows <- function(z) {
  notNanRow <- which(apply(dfClust[,stageNames], 1, function(x){sum(is.nan(x))}) == 0)
  return(notNanRow)
}

## Make heatmaps
makeHeatMapTableFromDeepToolsMatrix <- function(pM,
                                                xSmooth=5,
                                                ySmooth=20,
                                                midRange=20,
                                                binSize=10,
                                                winSize=5,
                                                name='ALL'){
  
  ## This is required to get ggplot to make nice looking heatmaps
  ## Horizontal AND vertical smoothing
  pM <- t(apply(pM,1,function(x){frollmean(x,xSmooth,align='center')}))
  pM <- apply(pM,2,function(x){frollmean(x,ySmooth,align='center')})
  
  # Take every nth line (reduce HM complexity)
  pM <- pM[seq(1,dim(pM)[1],ySmooth),]
  
  dfM <- as.data.frame(pM)
  names(dfM) <- 1:(dim(dfM)[2])
  
  ## Sort DF (max to min)
  dfM$score <- apply(dfM,1,FUN=function(x){nx=round(length(x)/2,0); mean(x[nx-midRange:nx+midRange])})
  dfS <- dfM[order(dfM$score),]
  
  ## Remove score column & add hs column
  dfS <- dfS[,1:(dim(dfS)[2]-1)]
  dfS$hs <- 1:dim(dfS)[1]
  
  mS <- melt(dfS,id.vars = 'hs')
  names(mS) <- c('hs','pos','score')
  
  ## Account for bin size in x-coordinates
  mS$pos <- as.numeric(mS$pos) - binSize/2
  mS$pos <- ((mS$pos-max(mS$pos)/2) / (max(mS$pos)/2) * winSize) 
  
  ## add label
  mS$label <- name
  
  ## Remove NAs
  mS <- mS[!is.na(mS$score),]
  
  return(mS)
}

plotFromDeepToolsMatrix <- function(inMatrix=NULL,
                                    inData=NULL,
                                    sortOrder=NULL,
                                    normalizeByIgG=FALSE,
                                    normByRow=TRUE,
                                    useAll=FALSE,
                                    matrixType=NULL,
                                    trimDate=TRUE,
                                    fSize=5,
                                    withHeatmap=FALSE,
                                    samplesToExclude=NULL){
  
  '%ni%'      <- Negate("%in%")
  
  if (is.null(matrixType)){return(0)}
  
  if (matrixType %ni% c('referencepoint','scaleregions')){
    print('Invalid matrixType [OPTS: referencepoint OR scaleregions]')
    return(FALSE)
  }
  
  print (inMatrix)
  m      <- fread(inMatrix,skip=1)
  m      <- m[,7:dim(m)[2]]
  
  sampleNames <- data.frame(name = sortOrder, 
                            initname = sortOrder, 
                            useMe = 1)
  
  mSize <- dim(m)[1]; mCol  <- dim(m)[2]
  
  nPerSet <- mCol/dim(sampleNames)[1]
  
  doneOne <- FALSE
  for (n in seq(nPerSet,mCol,by = nPerSet)){
    
    i        <- which(seq(nPerSet,mCol,by = nPerSet) == n)
    dStage <- sortOrder[i]
    
    sampleOK <- TRUE
    if (!useAll){
      if (sampleNames$useMe[(n/nPerSet)] == 0){
        sampleOK <- FALSE  
      }
    }
    
    if (!is.null(samplesToExclude)){
      if (sampleNames$name[(n/nPerSet)] %in% samplesToExclude){
        sampleOK <- FALSE  
      }
    }
    
    if (sampleOK){  
      thisData    <- m[,(n-(nPerSet-1)):n]
      if (withHeatmap){
        heatmapData <- makeHeatMapTableFromDeepToolsMatrix(thisData, binSize = 10, winSize = fSize, name= sampleNames$name[(n/nPerSet)])
      }else{
        heatmapData <- data.frame(x=NA,y=NA,label=NA)
      }
      
      nFrom    <- (1+mSize*((n/nPerSet)-1))
      nTo      <- mSize*((n/nPerSet))
      
      d <- t(apply(thisData,1,FUN=function(x){standardizeMNSD(x)}))
      
      ## Set NAs to 0
      d[is.na(d)] <- 0
      
      nFrom <- (n-(nPerSet-1))
      
      q33 <- as.integer(quantile(nFrom:n,c(5/12)))
      q67 <- as.integer(quantile(nFrom:n,c(7/12)))
      
      mMid <- as.matrix(m[,q33:(q67)])
      mL <- as.matrix(m[,nFrom:(q33)])
      mR <- as.matrix(m[,(q67):n])
      
      df <- data.frame(stage=dStage,
                       mid=rowSums(mMid),
                       flank=(rowSums(mL)+rowSums(mR))/2)
      
      ## Remove rows with no signal
      #d <- d[rowSums(d)>0,]
      
      #dS <- standardizeMNSD(d)
      #dS <- dS
      #cM <- colMeans(dS) - min(colMeans(dS))
      cM <- colMeans(d,na.rm=TRUE)
      cM <- cM - min(cM) + 0.1
      
      #fMean <- mean(cM[c(1:(nPerSet*0.10),(nPerSet*0.90):nPerSet)])
      #cM <- cM-fMean
      
      #sOrder$E[(n/nPerSet)] <- max(cM)  
      cD <- data.frame(pos    = seq(1,100,length.out = nPerSet),
                       FC     = cM,
                       normFC = cM/mean(cM[c(1:20,(nPerSet-20):nPerSet)]),
                       sample = sampleNames$name[(n/nPerSet)])
      
      cD$ok <- FALSE
      
      if (max(cD$FC) > 0.2){
        cD$ok <- TRUE  
      }
      
      if (!doneOne){
        allData  <- cD
        allHM    <- heatmapData
        distData <- df
        doneOne  <- TRUE
      }else{
        allData  <- rbind(allData,cD)
        allHM    <- rbind(allHM,heatmapData)
        distData <- rbind(distData,df)
      }
    }  
  }
  
  if (is.null(sortOrder)){
    ## Order by max enrichment
    sampleAgg <- aggregate(allData$FC,
                           by=list(allData$sample),
                           FUN=function(x){sum(x[x>quantile(x,.95)])})
    
    names(sampleAgg) <- c('sample','max')
    
    sampleOrder <- sampleAgg[order(sampleAgg$max),]
    
    
    allData$sample <- factor(allData$sample,
                             levels=sampleOrder$sample)
  }else{
    allData$sample <- factor(allData$sample,
                             levels=rev(sortOrder))
  }
  
  ## OK. Very specific. 
  ## NORMALIZE BY IgG
  if (normalizeByIgG){
    allData$normFactor <- allData$FC[allData$sample == 'IgG']
    allData$oldFC <- allData$FC
    allData$rawFC <- allData$oldFC-allData$normFactor
    
    allData$l2NormFactor <- allData$normFC[allData$sample == 'IgG']
    allData$FC           <- log2(allData$normFC/allData$l2NormFactor)
  }
  
  ## Exclude samples
  if (!is.null(samplesToExclude)){
    allData  <- allData[allData$sample %ni% samplesToExclude,]
    distData <- allData[allData$stage %ni% samplesToExclude,]
  }
  
  ## Midpoint is not really the midpoint ... it's 40% point
  ## ...... looks better
  #myMid <- (min(allData$FC)+max(allData$FC))*.2
  myMid <- 0
  
  if (matrixType == 'referencepoint'){
    
    g <- ggplot(allData,aes(x=pos,fill=FC,y=sample)) + 
      geom_tile(color=NA) + 
      scale_fill_gradient2(low='dodgerblue2',mid='white',high='red',midpoint = 0) +
      xlab('Distance to hotspot center (Kb)') + 
      ylab('')
    
    g <- g + scale_x_continuous(breaks=quantile(allData$pos,c(0,.5,1)),
                                labels=c(paste0('-',fSize),'0',paste0('+',fSize))) + 
      xlab('Position (Kbp)')
    
  }else{
    
    g <- ggplot(allData,aes(x=pos,fill=FC,y=sample)) + 
      geom_tile(color=NA) + 
      scale_fill_gradient2(low='dodgerblue2',mid='white',high='red',midpoint = 0) +
      xlab('Distance to hotspot center (Kb)') + 
      ylab('')
    
    g <- g + scale_x_continuous(breaks=quantile(allData$pos,c(0,.25,.75,1)),
                                labels=c(paste0('-',fSize),'TSS','TES',paste0('+',fSize))) + 
      xlab('Position (Kbp)')
    
  }  
  
  if (withHeatmap){
    
    allHM$label <- factor(allHM$label, levels=sortOrder)
    
    gHM <- ggplot(allHM,aes(x=pos,y=hs,fill=score)) + 
      geom_raster(na.rm=TRUE) + 
      scale_fill_gradientn(colors=c('grey20','white','firebrick'),
                           values=c(0,.4,1)) + 
      theme(legend.position="none",
            axis.line.y  = element_blank(),
            axis.ticks.y = element_blank(),
            axis.text.y  = element_blank(),
            strip.background=element_blank(),
            strip.text=element_text(face='bold',size=7)) + 
      coord_cartesian(xlim=c(-fSize*.9,fSize*.9),expand=FALSE) +
      ylab('') + 
      xlab('Dist to center (Kb)') + 
      scale_x_continuous(breaks=sort(unique(c(seq(0,fSize*.9,3)*-1,seq(0,fSize*.9,3))))) + 
      theme(plot.margin = unit(c(0,0,0,0),'cm')) + 
      guides(yLab = element_blank()) + facet_grid(~label)
    
    ## Average plots
    meanData <- aggregate(allHM$score,by=list(allHM$pos,allHM$label),FUN=mean)
    names(meanData) <- c('pos','name','score')
    
    ## Assure that ticks are pretty !!
    yN <- pretty(max(meanData$score))
    yMaxTick <- min(yN[yN != 0])
                    
    yN <- pretty(min(meanData$score))
    yMinTick <- max(yN[yN != 0])
    
    if (yMinTick >= (0 - 0.1*yMaxTick) & yMaxTick > 0){
      
      myMax <- max(c(yMaxTick,meanData$score))
      
      if (myMax < max(meanData$score) + 0.1*max(meanData$score)){
        yMaxVal  = max(meanData$score) + 0.1*max(meanData$score)
      }else{
        yMaxVal = myMax
      }
      
      yMinTick = 0;
      yMinVal  = 0 - 0.1*yMaxVal;
      
      yTicks <- c(0,yMaxTick)
    }else{
      if (yMinTick < 0 & yMaxTick < 0){
        myMin <- min(c(yMinTick,meanData$score))
        
        if (myMin > min(meanData$score) + 0.1*min(meanData$score)){
          yMinVal  = min(meanData$score) + 0.1*min(meanData$score)
        }else{
          yMinVal <- myMin
        }
        
        yMaxTick = 0;
        yMaxVal  = 0 + 0.1*yMinVal;
        yTicks <- c(yMinTick,0)
      }else{
          yMinVal = yMinTick - 0.1*max(meanData$score);
          yMaxVal = yMaxTick + 0.1*max(meanData$score);
          yTicks <- c(yMinTick,0,yMaxTick)
      }
    }
    
    print(paste0(yTicks))
    print(paste(yMinVal,yMaxVal,sep = " ... "))
    
    gMean <- ggplot(meanData,aes(x=pos,y=score)) + 
      geom_smooth(method='loess',span=.1,lwd=.2,color='dodgerblue2') + 
      theme(legend.position="none",
            axis.text.x=element_blank(),
            strip.background=element_blank(),
            strip.text=element_text(face='bold',size=7)) + 
      coord_cartesian(xlim=c(-fSize*.9,fSize*.9),
                      ylim=c(yMinVal,yMaxVal),
                      expand=FALSE) +
      ylab('Mean') + 
      xlab('') +
      scale_x_continuous(breaks=sort(unique(c(seq(0,fSize*.9,3)*-1,seq(0,fSize*.9,3))))) + 
      scale_y_continuous(breaks=yTicks) +
      theme(plot.margin = unit(c(0,0,0,0),'cm')) + 
      guides(yLab = element_blank()) + facet_grid(~name) 
      
    gHMandMean <- ggarrange(gMean,gHM,
                            ncol=1,nrow=2,
                            heights=c(2,8),
                            align='v')
      
  }else{
    gHM <- ggplot()
    gMean <- ggplot()
    gHMandMean <- ggplot()
    meanData <- NULL
  }

  return(list(fig=g,
              data=allData,
              order=levels(allData$sample),
              matrix=m,
              hmData = allHM, 
              figHM = gHM,
              figMean = gMean, 
              figBoth = gHMandMean,
              mData = meanData,
              dDist=distData))
}
