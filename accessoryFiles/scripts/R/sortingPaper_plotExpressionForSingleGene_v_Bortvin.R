#source('scripts/sortingPaper_functions.R')
#source('/share/code/R/RDCOstandardfunctions.R')

theme7point()

getData <- function(x){
  h3Exp   <- read.table('../../data/bigTable.withExp.tab',header=TRUE)
  ## Decide which dpp expression datapoints to use
  #expLabel <- c(12,14,16,18,20)
  #expNames <- paste0('d',expLabel)
  #stageLabel <- c("L","Z","E","P","D")
  #stageNames <- c("k4LE","k4ZY","k4EP","k4LP","k4DI")
  
  expLabel <- c('L','Z','P','D')
  expNames <- paste0('e',c('LE','ZY','PA','DI'))
  
  stageLabel <- c("L","Z","P","D")
  stageNames <- c("k4LE","k4ZY","k4PA","k4DI")
  
  ## Require at least 1.5 RPKM @ max timepoint
  #h3Exp$OK    <- apply(h3Exp[,expNames], 1, max) >= 0.5
  
  #h3Exp       <- h3Exp[h3Exp$OK,]
  
  h3Exp$d6  <- h3Exp$d6/sum(h3Exp$d6)*1000000
  h3Exp$d10 <- h3Exp$d10/sum(h3Exp$d10)*1000000
  h3Exp$d12 <- h3Exp$d12/sum(h3Exp$d12)*1000000
  h3Exp$d14 <- h3Exp$d14/sum(h3Exp$d14)*1000000
  h3Exp$d16 <- h3Exp$d16/sum(h3Exp$d16)*1000000
  h3Exp$d18 <- h3Exp$d18/sum(h3Exp$d18)*1000000
  h3Exp$d20 <- h3Exp$d20/sum(h3Exp$d20)*1000000
  #h3Exp$d38 <- h3Exp$d38/sum(h3Exp$d38)*1000000
  
  source('scripts/sortingPaper_getH3K4me3_and_H3K9ac_table.R')
  #source('scripts/sortingPaper_getH3K4me3table.R')
  ## h3 is from sortingPaper_getH3K4me3_and_H3K9ac.R (Above)
  ## take just TSS and just a few columns for the join with expression data
  '%ni%'        <- Negate("%in%")
  okH3$k4PA     <- (okH3$k4EP + okH3$k4LP)/2;
  justTSSs      <- okH3[okH3$type %ni% c('HS','HSTSS','HSKOTSS'),]
  h3Exp2merge   <- h3Exp[,c(1:3,43:55)]
  
  ## Join with expression data
  h4 <- join(justTSSs, h3Exp2merge, type='inner')
  
  ## Subset down to only columns of interest
  dTSS <- h4[,c("cs","from","to","type","name",
                stageNames,
                expNames)]
  return(dTSS)
}

plotGeneExp <- function(geneName = 'Hormad1',
                        dTSS = tssIn,
                        tall = TRUE){
  
  dMe <- dTSS[dTSS$name == geneName,]
  dMe$iso <- 1:length(dMe$cs)
  
  if (dim(dMe)[1]>0){
    m <- melt(dMe,id.vars = c('cs','from','to','type','name','iso'))
    
    m$dType[m$variable %in% c('eLE','eZY','ePA','eDI')] <- 'Expr'
    m$dType[m$variable %in% c('k4LE','k4ZY','k4PA','k4DI')] <- 'K4me3'
    
    m <- m[m$dType %in% c('Expr','K4me3'),]
    
    m$stage[m$variable %in% c('eLE','k4LE')] <- 'LE'
    m$stage[m$variable %in% c('eZY','k4ZY')] <- 'ZY'
    m$stage[m$variable %in% c('ePA','k4PA')] <- 'PA'
    m$stage[m$variable %in% c('eDI','k4DI')] <- 'DI'
    
    m$stage <- factor(m$stage,levels=c('LE','ZY','PA','DI'))
    
    g <- ggplot(m,aes(x=stage,y=as.numeric(value),group=iso,color=factor(iso))) + 
      facet_wrap(~dType,scales='free_y',ncol=1) + 
      geom_line() + geom_point() + 
      ggtitle(label=geneName)
  }else{
    print(paste0("No data for ",geneName))
    g<-NA
  }
  
  return(g)
}

plotGeneExpNF <- function(geneName = 'Hormad1',
                        dTSS = tssIn,
                        noXlab=FALSE,
                        noYlab=FALSE,
                        tall = TRUE){
  
  dMe <- dTSS[dTSS$name == geneName,]
  dMe$iso <- 1:length(dMe$cs)
  
  if (dim(dMe)[1]>0){
    m <- melt(dMe,id.vars = c('cs','from','to','type','name','iso'),
              measure.vars = c('eLE','eZY','ePA','eDI','k4LE','k4ZY','k4PA','k4DI'))
    
    m$dType[m$variable %in% c('eLE','eZY','ePA','eDI')] <- 'Expr'
    m$dType[m$variable %in% c('k4LE','k4ZY','k4PA','k4DI')] <- 'K4me3'
    
    m <- m[m$dType %in% c('Expr','K4me3'),]
    
    m$stage[m$variable %in% c('eLE','k4LE')] <- 'LE'
    m$stage[m$variable %in% c('eZY','k4ZY')] <- 'ZY'
    m$stage[m$variable %in% c('ePA','k4PA')] <- 'PA'
    m$stage[m$variable %in% c('eDI','k4DI')] <- 'DI'
    
    m$stage <- factor(m$stage,levels=c('LE','ZY','PA','DI'))
    
    g <- ggplot(m,aes(x=stage,
                      y=as.numeric(value),
                      group=dType,
                      color=dType)) + 
      geom_line() + geom_point() + 
      ggtitle(label=geneName)
    
  }else{
    print(paste0("No data for ",geneName))
    g<-NA
  }

  if (noXlab){
    g <- g + xlab('')
  }
  
  if (noYlab){
    g <- g + ylab('')
  }
  
  return(g)
}

plotGeneExpNFMulti <- function(geneNames = c('Dmc1','Mei1'),
                          dTSS = tssIn,
                          noXlab=FALSE,
                          noYlab=FALSE,
                          tall = TRUE){
  
  if (exists('dMe')){
    remove('dMe')
  }
  
  #print(geneNames)
  for (g in geneNames){
    #print(g)
    dM <- dTSS[dTSS$name == g,]
    dM$iso <- 1:length(dM$cs)
    
    #print(dim(dM))
    
    if (exists('dMe')){
      dMX <- dMe
      dMe <- rbind(dMX,dM)
    }else{
      dMe <- dM
    }
    #print(dim(dMe))
  }
  
  #print(dim(dMe))
  #print(dMe$name)
  
  if (dim(dMe)[1]>0){
    m <- melt(dMe,id.vars = c('cs','from','to','type','name','iso'),
              measure.vars = c('eLE','eZY','ePA','eDI','k4LE','k4ZY','k4PA','k4DI'))
    
    m$dType[m$variable %in% c('eLE','eZY','ePA','eDI')] <- 'RNA-Seq'
    m$dType[m$variable %in% c('k4LE','k4ZY','k4PA','k4DI')] <- 'H3K4me3'
    
    m$stage[m$variable %in% c('eLE','k4LE')] <- 'LE'
    m$stage[m$variable %in% c('eZY','k4ZY')] <- 'ZY'
    m$stage[m$variable %in% c('ePA','k4PA')] <- 'PA'
    m$stage[m$variable %in% c('eDI','k4DI')] <- 'DI'
    
    m$stage <- factor(m$stage,levels=c('LE','ZY','PA','DI'))
    
    if (tall){
    g <- ggplot(m,aes(x=stage,
                      y=as.numeric(value),
                      group=dType,
                      color=dType)) + 
      geom_line(lwd=.2) + geom_point(lwd=.2,size=.5) + 
      facet_wrap(~name,ncol=1) + 
      theme(strip.background=element_blank())
    }else{
      g <- ggplot(m,aes(x=stage,
                        y=as.numeric(value),
                        group=dType,
                        color=dType)) + 
        geom_line(lwd=.2) + geom_point(lwd=.2,size=.5) + 
        facet_wrap(~name,nrow=1) + 
        theme(strip.background=element_blank())
    }    
  }else{
    print(paste0("No data for ",geneName))
    g<-NA
  }
  
  if (noXlab){
    g <- g + xlab('')
  }
  
  if (noYlab){
    g <- g + ylab('')
  }
  
  return(g)
}

#myTSS <- getData()
#plotGeneExp('Hormad1',myTSS)
