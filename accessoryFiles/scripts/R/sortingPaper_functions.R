## Load default RDCO functions
source(paste0(Sys.getenv('SHARE'),'/code/R/RDCOstandardfunctions.R'))

## Set publication ready theme
theme7point()

## Define project-specific functions

standardize <- function(z) {
  return(standardizeMedMad(z))
}

nanRows <- function(z) {
  nanRow <- which(apply(z, 1, function(x){sum(is.nan(x))}) > 0)
  return(nanRow)
}

notNanRows <- function(z) {
  notNanRow <- which(apply(z, 1, function(x){sum(is.nan(x))}) == 0)
  return(notNanRow)
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
    print(paste0('Rep ',n,' ... tRem = ',round(nRem),' s'))
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

getStableTSS <- function(initMat, nLow = 0.15, nHi = 0.01, nFold = 1.5){
  nameFields <- c("PL1","LE1","LE2","ZY1","ZY2","ZY3","EP1","EP2","EP3","LP1","LP2","DI1","DI2")
  
  ## Remove cross reactors (HS @ TSS, TSS @ HS)
  h3Clean <- initMat[initMat$cleanType,]
  
  # # Remove below signal vals
  # h3Clean <- h3Clean[h3Clean$PL1 > 0,]
  # h3Clean <- h3Clean[h3Clean$LE1 > 0,]
  # h3Clean <- h3Clean[h3Clean$LE2 > 0,]
  # h3Clean <- h3Clean[h3Clean$ZY1 > 0,]
  # h3Clean <- h3Clean[h3Clean$ZY2 > 0,]
  # h3Clean <- h3Clean[h3Clean$ZY3 > 0,]
  # h3Clean <- h3Clean[h3Clean$EP1 > 0,]
  # h3Clean <- h3Clean[h3Clean$EP2 > 0,]
  # h3Clean <- h3Clean[h3Clean$EP3 > 0,]
  # h3Clean <- h3Clean[h3Clean$LP1 > 0,]
  # h3Clean <- h3Clean[h3Clean$LP2 > 0,]
  # h3Clean <- h3Clean[h3Clean$DI1 > 0,]
  # h3Clean <- h3Clean[h3Clean$DI2 > 0,]

  ## Split to TSS and HS
  hTSS <- h3Clean[h3Clean$type %in% c('TSS','TSSKO'),]

  ## Establish boundaries
  pllo   <- quantile(hTSS$PL1,nLow,na.rm=TRUE)
  le1lo  <- quantile(hTSS$LE1,nLow,na.rm=TRUE)
  le2lo  <- quantile(hTSS$LE2,nLow,na.rm=TRUE)
  zy1lo  <- quantile(hTSS$ZY1,nLow,na.rm=TRUE)
  zy2lo  <- quantile(hTSS$ZY2,nLow,na.rm=TRUE)
  zy3lo  <- quantile(hTSS$ZY3,nLow,na.rm=TRUE)
  ep1lo  <- quantile(hTSS$EP1,nLow,na.rm=TRUE)
  ep2lo  <- quantile(hTSS$EP2,nLow,na.rm=TRUE)
  ep3lo  <- quantile(hTSS$EP3,nLow,na.rm=TRUE)
  lp1lo  <- quantile(hTSS$LP1,nLow,na.rm=TRUE)
  lp2lo  <- quantile(hTSS$LP2,nLow,na.rm=TRUE)
  di1lo  <- quantile(hTSS$DI1,nLow,na.rm=TRUE)
  di2lo  <- quantile(hTSS$DI2,nLow,na.rm=TRUE)
  
  plhi   <- quantile(hTSS$PL1,1-nHi,na.rm=TRUE)
  le1hi  <- quantile(hTSS$LE1,1-nHi,na.rm=TRUE)
  le2hi  <- quantile(hTSS$LE2,1-nHi,na.rm=TRUE)
  zy1hi  <- quantile(hTSS$ZY1,1-nHi,na.rm=TRUE)
  zy2hi  <- quantile(hTSS$ZY2,1-nHi,na.rm=TRUE)
  zy3hi  <- quantile(hTSS$ZY3,1-nHi,na.rm=TRUE)
  ep1hi  <- quantile(hTSS$EP1,1-nHi,na.rm=TRUE)
  ep2hi  <- quantile(hTSS$EP2,1-nHi,na.rm=TRUE)
  ep3hi  <- quantile(hTSS$EP3,1-nHi,na.rm=TRUE)
  lp1hi  <- quantile(hTSS$LP1,1-nHi,na.rm=TRUE)
  lp2hi  <- quantile(hTSS$LP2,1-nHi,na.rm=TRUE)
  di1hi  <- quantile(hTSS$DI1,1-nHi,na.rm=TRUE)
  di2hi  <- quantile(hTSS$DI2,1-nHi,na.rm=TRUE)
  
  ## Get TSSs within boundaries
  hTSS$PL1ok  <- TRUE
  hTSS$LE1ok  <- TRUE
  hTSS$LE1ok  <- TRUE
  hTSS$ZY1ok  <- TRUE
  hTSS$ZY2ok  <- TRUE
  hTSS$ZY3ok  <- TRUE
  hTSS$EP1ok  <- TRUE
  hTSS$EP3ok  <- TRUE
  hTSS$EP3ok  <- TRUE
  hTSS$LP1ok  <- TRUE
  hTSS$LP2ok  <- TRUE
  hTSS$DI1ok  <- TRUE
  hTSS$DI2ok  <- TRUE
  
  hTSS$ok  <- TRUE
  
  hTSS$ok[ (hTSS$PL1 >= plhi)  | (hTSS$PL1 <= pllo)]   <- FALSE
  hTSS$ok[ (hTSS$LE1 >= le1hi) | (hTSS$LE1 <= le1lo)]  <- FALSE
  hTSS$ok[ (hTSS$LE2 >= le2hi) | (hTSS$LE2 <= le2lo)]  <- FALSE
  hTSS$ok[ (hTSS$ZY1 >= zy1hi) | (hTSS$ZY1 <= zy1lo)]  <- FALSE
  hTSS$ok[ (hTSS$ZY2 >= zy2hi) | (hTSS$ZY2 <= zy2lo)]  <- FALSE
  hTSS$ok[ (hTSS$ZY3 >= zy3hi) | (hTSS$ZY3 <= zy3lo)]  <- FALSE
  hTSS$ok[ (hTSS$EP1 >= ep1hi) | (hTSS$EP1 <= ep1lo)]  <- FALSE
  hTSS$ok[ (hTSS$EP2 >= ep2hi) | (hTSS$EP1 <= ep2lo)]  <- FALSE
  hTSS$ok[ (hTSS$EP3 >= ep3hi) | (hTSS$EP1 <= ep3lo)]  <- FALSE
  hTSS$ok[ (hTSS$LP1 >= lp1hi) | (hTSS$LP1 <= lp1lo)]  <- FALSE
  hTSS$ok[ (hTSS$LP2 >= lp2hi) | (hTSS$LP2 <= lp2lo)]  <- FALSE
  hTSS$ok[ (hTSS$DI1 >= di1hi) | (hTSS$DI1 <= di1lo)]  <- FALSE
  hTSS$ok[ (hTSS$DI2 >= di2hi) | (hTSS$DI2 <= di2lo)]  <- FALSE
  
  ## Make "final" set of TSSs to study
  okTSS <- hTSS[hTSS$ok, ]
  
  okTSS$allOK <- TRUE
  
  for (n1 in nameFields){
    for (n2 in nameFields){
      ## Sum of expression
      okTSS[,paste0('tot_',n1,'_v_',n2)] <- (okTSS[,n1]+okTSS[,n2])/2
      
      ## Log2 difference - median(log2 diff)
      okTSS[,paste0('l2_',n1,'_v_',n2)]  <- log2(okTSS[,n1]/okTSS[,n2])  - log2(median(okTSS[,n1]/okTSS[,n2]))
      
      ## replace NaNs with 0s
      okTSS[is.nan(okTSS[,paste0('l2_',n1,'_v_',n2)]),paste0('l2_',n1,'_v_',n2)] <- 0
      
      ## Mark TSS with log2(delta) > log2(nFold)
      okTSS$allOK[abs(okTSS[,paste0('l2_',n1,'_v_',n2)]) >= log2(nFold)]  <- FALSE
    }
  }
  
  return(okTSS)
}

getStableTSS_new <- function(initMat, nLow = 0.15, nHi = 0.01, nFold = 1.2,modType='h3k4me3'){
  
  if (modType == 'h3k4me3'){
    initMat$LE <- initMat$k4LE
    initMat$ZY <- initMat$k4ZY
    initMat$EP <- initMat$k4EP
    initMat$LP <- initMat$k4LP
    initMat$DI <- initMat$k4DI
  }else{
    initMat$LE <- initMat$k9LE
    initMat$ZY <- initMat$k9ZY
    initMat$EP <- initMat$k9EP
    initMat$LP <- initMat$k9LP
    initMat$DI <- initMat$k9DI
  }
  
  nameFields <- c("LE","ZY","EP","LP","DI")
  
  ## Remove cross reactors (HS @ TSS, TSS @ HS)
  h3Clean <- initMat[initMat$cleanType,]
  
  # # Remove below signal vals
  for (n in nameFields){
    h3Clean <- h3Clean[h3Clean[,n] > 0,]  
  }

  ## Split to TSS and HS
  hTSS <- h3Clean[h3Clean$type %in% c('TSS','TSSKO'),]
  
  ## Establish boundaries
  LElo  <- quantile(hTSS$LE,nLow,na.rm=TRUE)
  ZYlo  <- quantile(hTSS$ZY,nLow,na.rm=TRUE)
  EPlo  <- quantile(hTSS$EP,nLow,na.rm=TRUE)
  LPlo  <- quantile(hTSS$LP,nLow,na.rm=TRUE)
  DIlo  <- quantile(hTSS$DI,nLow,na.rm=TRUE)
  
  LEhi  <- quantile(hTSS$LE,1-nHi,na.rm=TRUE)
  ZYhi  <- quantile(hTSS$ZY,1-nHi,na.rm=TRUE)
  EPhi  <- quantile(hTSS$EP,1-nHi,na.rm=TRUE)
  LPhi  <- quantile(hTSS$LP,1-nHi,na.rm=TRUE)
  DIhi  <- quantile(hTSS$DI,1-nHi,na.rm=TRUE)
  
  ## Get TSSs within boundaries
  hTSS$LEok   <- TRUE
  hTSS$ZYok   <- TRUE
  hTSS$EPok   <- TRUE
  hTSS$LPok   <- TRUE
  hTSS$DIok   <- TRUE
  
  hTSS$ok  <- TRUE
  
  hTSS$ok[ (hTSS$LE >= LEhi)  | (hTSS$LE <= LElo)]   <- FALSE
  hTSS$ok[ (hTSS$ZY >= ZYhi)  | (hTSS$ZY <= ZYlo)]   <- FALSE
  hTSS$ok[ (hTSS$EP >= EPhi)  | (hTSS$EP <= EPlo)]   <- FALSE
  hTSS$ok[ (hTSS$LP >= LPhi)  | (hTSS$LP <= LPlo)]   <- FALSE
  hTSS$ok[ (hTSS$DI >= DIhi)  | (hTSS$DI <= DIlo)]   <- FALSE
  
  ## Make "final" set of TSSs to study
  okTSS <- hTSS[hTSS$ok, ]
  
  okTSS$allOK <- TRUE
  
  for (n1 in nameFields){
    for (n2 in nameFields){
      ## Sum of expression
      okTSS[,paste0('tot_',n1,'_v_',n2)] <- (okTSS[,n1]+okTSS[,n2])/2
      
      ## Log2 difference - median(log2 diff)
      okTSS[,paste0('l2_',n1,'_v_',n2)]  <- log2(okTSS[,n1]/okTSS[,n2])  - log2(median(okTSS[,n1]/okTSS[,n2]))
      
      ## replace NaNs with 0s
      okTSS[is.nan(okTSS[,paste0('l2_',n1,'_v_',n2)]),paste0('l2_',n1,'_v_',n2)] <- 0
      
      ## Mark TSS with log2(delta) > log2(nFold)
      okTSS$allOK[abs(okTSS[,paste0('l2_',n1,'_v_',n2)]) >= log2(nFold)]  <- FALSE
    }
  }
  
  return(okTSS)
}

correctH3K4 <- function(initMat, mTSS){
  nameFields <- c("PL1","LE1","LE2","ZY1","ZY2","ZY3","EP1","EP2","EP3","LP1","LP2","DI1","DI2")
  nameFields <- c("LE1","LE2","ZY1","ZY2","ZY3","EP1","EP2","EP3","LP1","LP2","DI1","DI2")
 
  cFactor <- data.frame(sample     = nameFields,
                        corrFactor = 0)
  
  for (n1 in nameFields){
    correctionFactor <- median(mTSS[mTSS$allOK,n1])
    initMat[,n1] <- initMat[,n1]/correctionFactor
    cFactor[cFactor$sample == n1,2] <- correctionFactor
  }
  
  return(list(h3Mat=initMat,corrFactors=cFactor))
}

correctH3K9 <- function(initMat, mTSS){
  nameFields <- c("k9LE","k9ZY","k9EP","k9LP","k9DI")
  
  cFactor <- data.frame(sample     = nameFields,
                        corrFactor = 0)
  
  for (n1 in nameFields){
    correctionFactor <- median(mTSS[mTSS$allOK,n1])
    initMat[,n1] <- initMat[,n1]/correctionFactor
    cFactor[cFactor$sample == n1,2] <- correctionFactor
  }
  
  return(list(h3Mat=initMat,corrFactors=cFactor))
}

correctH3K4andH3K9 <- function(initMat, mK4TSS, mK9TSS){
  nameFields <- c("k9LE","k9ZY","k9EP","k9LP","k9DI")
  
  cFactorK9 <- data.frame(sample     = nameFields,
                        corrFactor = 0)
  
  for (n1 in nameFields){
    n1TPM                               <- paste0(n1,'TPM')
    correctionFactor                    <- median(mK9TSS[mK9TSS$allOK,n1])
    initMat[,n1]                        <- initMat[,n1]/correctionFactor
    initMat[,n1TPM]                     <- toTPM(initMat[,n1])
    cFactorK9[cFactorK9$sample == n1,2] <- correctionFactor
  }
  
  nameFields <- c("k4PL","k4LE","k4ZY","k4EP","k4LP","k4DI")
  
  cFactorK4 <- data.frame(sample     = nameFields,
                        corrFactor = 0)
  
  for (n1 in nameFields){
    n1TPM                               <- paste0(n1,'TPM')
    correctionFactor                    <- median(mK4TSS[mK4TSS$allOK,n1])
    initMat[,n1]                        <- initMat[,n1]/correctionFactor
    initMat[,n1TPM]                     <- toTPM(initMat[,n1])
    cFactorK4[cFactorK4$sample == n1,2] <- correctionFactor
  }
  
  cFactor <- rbind(cFactorK4,cFactorK9)
  
  return(list(h3Mat=initMat,corrFactors=cFactor))
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

getStageColors <- function(nAlpha=1,
                           darkLvl=1){
  
  cLE <- alpha(darkenColor('#c87137ff',darkLvl),nAlpha)
  cZY <- alpha(darkenColor('#ffb380ff',darkLvl),nAlpha)
  cEP <- alpha(darkenColor('#ffcc00ff',darkLvl),nAlpha)
  cLP <- alpha(darkenColor('#ffaaeeff',darkLvl),nAlpha)
  cDI <- alpha(darkenColor('#619cffff',darkLvl),nAlpha)
  
  return(c(cLE,cZY,cEP,cLP,cDI))
}


plotFromDeepToolsMatrix <- function(inMatrix=NULL,
                                    inData=NULL,
                                    sortOrder=NULL,
                                    normalizeByIgG=FALSE,
                                    samples='../../data/histoneCodeAtHS/alDT.names.tab',
                                    originalOrder=FALSE,
                                    useAll=FALSE,
                                    matrixType=NULL,
                                    trimDate=TRUE,
                                    fSize=5){
  
  '%ni%'      <- Negate("%in%")
  
  if (is.null(matrixType)){return(0)}
  
  if (matrixType %ni% c('referencepoint','scaleregions')){
    print('Invalid matrixType [OPTS: referencepoint OR scaleregions]')
    return(FALSE)
  }
  
  if (is.null(inData)){
    m      <- fread(inMatrix,skip=1)
    m      <- m[,7:dim(m)[2]]
  }else{
    m      <- inData
  }
  
  sampleNames <- read.table(samples)
  names(sampleNames) <- c('initname','useMe')
  sampleNames$name <- gsub(pattern = '_\\d{6}',
                           x = sampleNames$initname, 
                           replacement = '', 
                           perl=TRUE)
  
  mSize <- dim(m)[1]
  mCol  <- dim(m)[2]
  
  nPerSet <- mCol/dim(sampleNames)[1]
  
  doneOne <- FALSE
  
  #for (n in seq(nPerSet,mCol,by = nPerSet)){
  for (n in seq(nPerSet,mCol,by = nPerSet)){
    
    sampleOK <- TRUE
    if (!useAll){
      if (sampleNames$useMe[(n/nPerSet)] == 0){
        sampleOK <- FALSE  
      }
    }
    
    if (sampleOK){  
      thisData <- m[,(n-(nPerSet-1)):n]
      nFrom    <- (1+mSize*((n/nPerSet)-1))
      nTo      <- mSize*((n/nPerSet))
      
      d <- thisData
      
      ## Set NAs to 0
      d[is.na(d)] <- 0
      
      ## Remove rows with no signal
      d <- d[rowSums(d)>0,]
      
      dS <- standardizeMNSD(d)
      dS <- dS
      cM <- colMeans(dS) - min(colMeans(dS))
      
      fMean <- mean(cM[c(1:(nPerSet*0.10),(nPerSet*0.90):nPerSet)])
      cM <- cM-fMean
      
      #sOrder$E[(n/nPerSet)] <- max(cM)  
      cD <- data.frame(pos = seq(1,100,length.out = nPerSet),
                       FC = cM,
                       sample = sampleNames$name[(n/nPerSet)])
      
      cD$ok <- FALSE
      
      if (max(cD$FC) > 0.2){
        cD$ok <- TRUE  
      }
      
      if (!doneOne){
        allData <- cD
        doneOne <- TRUE
      }else{
        allData <- rbind(allData,cD)
      }
    }  
  }
  
  if (!originalOrder){
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
                               levels=sortOrder)
    }
  }
  
  ## Midpoint is not really the midpoint ... it's 40% point
  ## ...... looks better
  myMid <- (min(allData$FC)+max(allData$FC))*.4
  
  ## OK. Very specific. 
  ## NORMALIZE BY IgG
  if (normalizeByIgG){
    allData$normFactor <- allData$FC[allData$sample == 'IgG']
    
    allData$oldFC <- allData$FC
    allData$FC <- allData$FC-allData$normFactor
  }
  
  if (matrixType == 'referencepoint'){
    
    g <- ggplot(allData,aes(x=pos,fill=FC,y=sample)) + 
      geom_tile(color=NA) + 
      scale_fill_gradient2(low='grey80',
                           mid='firebrick',
                           high='yellow',
                           midpoint=myMid)  + 
      xlab('Distance to hotspot center (Kb)') + 
      ylab('')
    
    g <- g + scale_x_continuous(breaks=quantile(allData$pos,c(0,.5,1)),
                                labels=c(paste0('-',fSize),'0',paste0('+',fSize))) + 
      xlab('Position (Kbp)')
    
  }else{
    
    g <- ggplot(allData,aes(x=pos,fill=FC,y=sample)) + 
      geom_tile(color=NA) + 
      scale_fill_gradient2(low='grey80',mid='firebrick',high='yellow',midpoint=myMid)  + 
      xlab('Distance to hotspot center (Kb)') + 
      ylab('')
    
    g <- g + scale_x_continuous(breaks=quantile(allData$pos,c(0,.25,.75,1)),
                                labels=c(paste0('-',fSize),'TSS','TES',paste0('+',fSize))) + 
      xlab('Position (Kbp)')
    
  }  
  
  return(list(fig=g,
              data=allData,
              order=levels(allData$sample)))
}


plotHistoneHeatMaps <- function(justTop2k = TRUE,
                                justHS=FALSE,
                                SO=NULL,
                                usePaperSO=FALSE,
                                normByIgG=TRUE){
  theme7point() ## SET THEME
  
  if (justTop2k){
    top2k <- 'top2k_'
  }else{
    top2k <- ''
  }
  
  if (usePaperSO){
    SO <- rev(c('H3K4me3',
                'H3K9ac',
                'H3K4me2',
                'H3K36me3',
                'H4ac5',
                'H3K4me1',
                'H3K27ac',
                'H4K8ac',
                'H4K12ac',
                'H3K27me3',
                'H4K20me3',
                'H3K79me1',
                'H3K4ac',
                'H3K79me3',
                'H3K27me1',
                'H3K9me2',
                'H3K9me3',
                'H3',
                'IgG',
                'input'))
  }
  
  pScaleHS  <- plotFromDeepToolsMatrix(inMatrix = paste0('../../data/histoneCodeAtHS/deeptools/',top2k,'GL_histMods_July_2018_5k_plot3.v_B6oneMotif_HS.matrix'),
                                       matrixType = 'referencepoint',
                                       normalizeByIgG = normByIgG,
                                       sortOrder = SO)
  
  if (justHS){
    return(list(gHS=pScaleHS))
  }else{
    pScaleGene <- plotFromDeepToolsMatrix(inMatrix=paste0('../../data/histoneCodeAtHS/deeptools/',top2k,'GL_histMods_July_2018_5k_plot6.v_ScaledTranscripts.matrix'),
                                          sortOrder=pScaleHS$order,
                                          matrixType = 'scaleregions',
                                          normalizeByIgG = normByIgG)
    
    pScaleLINE <- plotFromDeepToolsMatrix(inMatrix=paste0('../../data/histoneCodeAtHS/deeptools/',top2k,'GL_histMods_July_2018_5k_plot7.v_ScaledL1LINEs.matrix'),
                                          sortOrder=pScaleHS$order,
                                          matrixType = 'scaleregions',
                                          normalizeByIgG = normByIgG)
    
    pScaleTSS <- plotFromDeepToolsMatrix(inMatrix = paste0('../../data/histoneCodeAtHS/deeptools/',top2k,'GL_histMods_July_2018_5k_plot4.v_TSS.matrix'),
                                         sortOrder=pScaleHS$order,
                                         matrixType = 'referencepoint',
                                         normalizeByIgG = normByIgG)
    
    pScaleTES <- plotFromDeepToolsMatrix(inMatrix = paste0('../../data/histoneCodeAtHS/deeptools/',top2k,'GL_histMods_July_2018_5k_plot5.v_TES.matrix'),
                                         sortOrder=pScaleHS$order,
                                         matrixType = 'referencepoint',
                                         normalizeByIgG = normByIgG)
    
    pScaleSINE <- plotFromDeepToolsMatrix(inMatrix = paste0('../../data/histoneCodeAtHS/deeptools/',top2k,'GL_histMods_July_2018_5k_plot12.v_SINE_Center.matrix'),
                                          sortOrder=pScaleHS$order,
                                          matrixType = 'referencepoint',
                                          normalizeByIgG = normByIgG)
    
    pScaleEnhancer <- plotFromDeepToolsMatrix(inMatrix = paste0('../../data/histoneCodeAtHS/deeptools/',top2k,'GL_histMods_July_2018_5k_plot11.v_Enhancers_RefSeq.matrix'),
                                              sortOrder=pScaleHS$order,
                                              matrixType = 'referencepoint',
                                              normalizeByIgG = normByIgG)
    
    noLeg <- theme(legend.position='none')
    noTxt <- theme(legend.position='none',axis.text.y=element_blank())
    
    gAll <- ggarrange(pScaleHS$fig +noLeg + ggtitle('HS'),
                      pScaleTSS$fig+noLeg+noTxt + ggtitle('TSS'),
                      pScaleTES$fig+noLeg+noTxt + ggtitle('TES'),
                      pScaleGene$fig+noLeg+noTxt + ggtitle('Genes'),
                      pScaleLINE$fig+noLeg+noTxt + ggtitle('L1 LINEs'),
                      pScaleSINE$fig+noLeg+noTxt + ggtitle('SINEs'),
                      pScaleEnhancer$fig+noLeg+noTxt + ggtitle('Enhancer'),
                      ncol=7,nrow=1,
                      widths=c(2,1,1,1,1,1,1))
    
    graphics.off()
    if (justTop2k){
      png(getIMGname('sortingPaper_allHistoneMarks_v_Genomic_Features_top2k','PNG'),width=10,height=4,units='in',res=300)
    }else{
      png(getIMGname('sortingPaper_allHistoneMarks_v_Genomic_Features_ALL','PNG'),width=10,height=4,units='in',res=300)
    }
    gAll
    dev.off()
    
    return(list(gHS =pScaleHS$fig,
                gTSS =pScaleTSS$fig,
                gTES =pScaleTES$fig,
                gGene =pScaleGene$fig,
                gLINE =pScaleLINE$fig,
                gSINE =pScaleSINE$fig,
                gEnh =pScaleEnhancer$fig,
                gall = gAll))
  }
  
}


drawSliceFig2B <- function(sliceNum=8,
                           featureSize=NULL,
                           corrFactors=NULL){
  
  library(png)
  
  s <- sliceNum
  slice <- sliceNum
    
  grHS  <- read.table(paste0('../../data/H3K4me3_x5/bam/dataSlice1/slice',slice,'.HSB6.bed'))
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
  
  grTSS  <- read.table(paste0('../../data/H3K4me3_x5/bam/dataSlice1/slice',slice,'.TSSrefseq.bed'))
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
  
  #colz <- brewer.pal(6, "RdYlGn")
  colfunc<-colorRampPalette(c("gold2","darkorange","firebrick4","grey30","forestgreen","springgreen"))
  colz <- colfunc(6)
  
  colz<-getStageColors()
  
  #zNorm <- load('../../data/H3K4me3_x5/Rtable/allH3K4me3_correctionFactors.R')
  zNorm <- load('../../data/H3K4me3_x5/Rtable/k4k9correctionFactors.R')
  
  #aPL <- read.table('../../data/H3K4me3_x5/bam/dataSlice1/slice1_PL_rep1_H3K4me3.mm10.bedgrap')
  aLE <- read.table(paste0('../../data/H3K4me3_x5/bam/dataSlice1/slice',slice,'_LE_rep2_H3K4me3.mm10.DZ.bedgraph'))
  aZY <- read.table(paste0('../../data/H3K4me3_x5/bam/dataSlice1/slice',slice,'_ZY_rep3_H3K4me3.mm10.DZ.bedgraph'))
  aEP <- read.table(paste0('../../data/H3K4me3_x5/bam/dataSlice1/slice',slice,'_EP_rep1_H3K4me3.mm10.DZ.bedgraph'))
  aLP <- read.table(paste0('../../data/H3K4me3_x5/bam/dataSlice1/slice',slice,'_LP_rep1_H3K4me3.mm10.DZ.bedgraph'))
  aDI <- read.table(paste0('../../data/H3K4me3_x5/bam/dataSlice1/slice',slice,'_DI_rep1_H3K4me3.mm10.DZ.bedgraph'))
  
  cf <- H3K4me3_GLcorrectionfactors
  ccLE <- cf[cf$sample == 'k4LE',2]/1000
  ccZY <- cf[cf$sample == 'k4ZY',2]/1000
  ccEP <- cf[cf$sample == 'k4EP',2]/1000
  ccLP <- cf[cf$sample == 'k4LP',2]/1000
  ccDI <- cf[cf$sample == 'k4DI',2]/1000
  
  aLE$V3 <- aLE$V3/ccLE
  aZY$V3 <- aZY$V3/ccZY
  aEP$V3 <- aEP$V3/ccEP
  aLP$V3 <- aLP$V3/ccLP
  aDI$V3 <- aDI$V3/ccDI
  
  aLE$stage <- 'LE'
  aZY$stage <- 'ZY'
  aEP$stage <- 'EP'
  aLP$stage <- 'LP'
  aDI$stage <- 'DI'
  
  cAll <- rbind(aLE,aZY,aEP,aLP,aDI)
  
  names(cAll) <- c('cs','pos','cover','stage')
  
  cAll$stage <- factor(cAll$stage,levels=c('LE','ZY','EP','LP','DI'))
  gInit <- ggplot(cAll,
                  aes(x=pos/1000000,
                      y=cover,
                      fill=stage,
                      color=stage))  
  
  for (hs in 1:dim(grHS)[1]){
    gHSfill <-annotate(geom='rect',
                       xmin=grHS$minX/1000000,
                       xmax=grHS$maxX/1000000,
                       ymin=-Inf,
                       ymax=Inf,
                       fill=alpha('black',.06),
                       color=NA)
  }
  
  for (hs in 1:dim(grTSS)[1]){
    gTSSfill <- annotate(geom='rect',
                         xmin=grTSS$minX/1000000,
                         xmax=grTSS$maxX/1000000,
                         ymin=-Inf,
                         ymax=Inf,
                         fill=alpha('red',.06),
                         color=NA)
  }
  
  diff <- (max(cAll$pos)-min(cAll$pos))*0.05
  xMin <- (min(cAll$pos)+diff)/1000000
  xMax <- (max(cAll$pos)-diff)/1000000
  
  getMaxVal <- function(x){
    floor(x*10^abs(floor(log10(x))))/10^abs(floor(log10(x)))
  }
  
  #yMax <- getMaxVal(max(cAll$cover))
  yMax <- round(max(cAll$cover),1)
  
  gData <- gInit + 
    gHSfill + gTSSfill + 
    geom_polygon(lwd=0.2) + 
    facet_grid(stage~.) + 
    coord_cartesian(xlim=c(xMin,xMax)) +
    scale_color_manual(values=getStageColors(darkLvl=1.7)) + 
    scale_fill_manual(values=getStageColors()) + 
    xlab('') + 
    ylab('H3K4me3 ChIP-Seq Coverage') + 
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
  
  gAnnotation <- ggplot(cAnnotation,
                        aes(xmin=minX/1000000,
                            xmax=maxX/1000000,
                            ymin=minY,
                            ymax=maxY,
                            fill=type,
                            color=type)) + 
    geom_rect(lwd=0.2) + 
    facet_grid(type~.) + 
    coord_cartesian(xlim=c(xMin,xMax),ylim=c(0.95,2.05)) +
    scale_color_manual(values=c(alpha(darkenColor('grey30',1),1),
                                alpha(darkenColor('red',1))),1) + 
    scale_fill_manual(values=c(alpha(darkenColor('grey30',1.3),1),
                               alpha(darkenColor('red',1.3),1))) + 
    xlab(paste0('Position on ',cAll$cs[1],' (Mbp)')) + 
    ylab('') + 
    theme(legend.position='none',
          axis.text.y  = element_blank(),        
          axis.ticks.y = element_blank(),
          axis.line.y = element_blank(),
          strip.background = element_blank(),
          strip.text = element_blank(),
          plot.margin = unit(c(0,0,0,0),'cm')) + 
    geom_text(aes(x, y, label=lab),
              data=data.frame(x=Inf, y=1.5, 
                              minX=Inf,maxX=Inf,
                              minY=Inf,maxY=Inf,
                              lab=c("Hotspots","TSS"),
                              type=c("Hotspots","TSS")),
              vjust=0.5, hjust=1, size=7*5/14)
  
  gCoverage <-
    ggarrange(gData,gAnnotation, 
                         align = 'v',
                         heights=c(10,2),
                         ncol=1,nrow=2)

  return(list(gFig=gCoverage,
              gTracks=gData,
              gAnno=gAnnotation)) 
}

plotPCAClusters <- function(myR, 
                            myI,
                            pc2use = c(1:4),
                            randomSeed=FALSE,
                            nClusters = 6,
                            histmod2use=NULL,
                            col=NULL){
  
  ## Clustering Plots
  ## make Cluster Plots
  myR$name <- factor(myR$name,levels=unique(myR$name))
  
  if (randomSeed){
    myCenters <- nClusters
  }else{
    nClusters <- length(levels(myR$name))
    cCenters <- vector(length = nClusters)
    for (nV in 1:length(cCenters)){
      cCenters[nV] <- min(which(myR$name==levels(myR$name)[nV]))
    }

    if (is.null(histmod2use)){
      myCenters <- myI$coord[cCenters,pc2use]
    }else{
      for (h in histmod2use){
        myR[[h]] <- standardizeMNSD(myR[[h]])
      }
      myCenters <- myR[cCenters,histmod2use]
    }
  }  
  
  if (is.null(histmod2use)){
    k      <- kmeans(myI$coord[,pc2use],
                    centers = myCenters,
                    iter.max = 40000)
  }else{
    k      <- kmeans(myR[,histmod2use],
                     centers = myCenters,
                     iter.max = 40000)
  }
  
  myR$k <- k$cluster
  
  ## Add N values to cluster names
  tN        <- as.data.frame(table(myR$k))
  names(tN) <- c('k','Freq')
  tN$kName <- paste0(tN$k,'; N=',format(x = tN$Freq,big.mark = ','))
  rK <- plyr::join(myR,tN,by='k')
  
  ## Group by cluster then by type (for return)
  rK$newName <- 'TSS'
  rK$newName[rK$name %in% c('Hot HS','HS')] <- 'HS' 
  myTableX <- table(rK[,c('kName','newName')])
  xTable     <- as.data.frame(prop.table(myTableX,1)*100)

  d           <- cbind(xTable,as.data.frame(myTableX))
  dAll        <- cbind(d[1:nClusters,],
                       d[(nClusters+1):(nClusters+nClusters),])[,c(1,3,6,9,12)]
  
  names(dAll) <- c('k','fHS','nHS','fTSS','nTSS')

  if (is.null(histmod2use)){
    dAll$type <- paste0('PCs: ',paste0(pc2use,'; ',collapse = "") )
    dAll$PCA  <- TRUE
    dAll$nProps <- length(pc2use)
  }else{
    dAll$type <- paste0('HM: ',paste0(histmod2use,'; ',collapse = "") )
    dAll$PCA  <- FALSE
    dAll$nProps <- length(histmod2use)
  }
  
  ## Group by cluster then by type
  myTable1 <- table(rK[,c('kName','name')])
  kTable     <- as.data.frame(prop.table(myTable1,1)*100)
  
  ## Group by type then by cluster
  myTable2 <- table(rK[,c('name','kName')])
  kTableCol  <- as.data.frame(prop.table(myTable2,1)*100)
  kTableColN <- as.data.frame(myTable2)
  
  ## PLOTS
  kTable$kName <- factor(kTable$kName,levels = rev(unique(kTable$kName)))
  myTheme <- theme(legend.title=element_blank())
  
  if (is.null(col)){
    myFill <- scale_shape_manual()
  }else{
    myFill <- scale_fill_manual(values=col)
  }
  
  gCKB <- ggplot(kTable,
                 aes(x=kName,
                     fill=name,
                     y=Freq)) + 
    geom_bar(stat='identity',width=.75) +  
    xlab('Cluster') + 
    ylab('Regions in cluster (%)') + 
    myFill + 
    scale_y_continuous(breaks=c(seq(0,100,20))) +
    geom_hline(yintercept=seq(20,80,20),lty='dotted',lwd=.2) +
    geom_hline(yintercept=seq(10,100,10),lty='dotted',lwd=.1) +
    coord_flip(ylim = c(0,100),
               xlim=c(0.5,nClusters+0.5),
               expand = FALSE) + 
    myTheme + 
    theme(panel.grid.major.y = element_line(linetype='dashed',
                                            size=.3),
          legend.position='bottom',
          legend.key.size=unit(0.2,'cm'),
          legend.text=element_text(size=7),
          axis.text.y = element_blank()) + 
    geom_label(aes(label=kName),
               y=50,fill='white',
               size=7*5/14)
  
  gCC <- ggplot(kTableCol, aes(x=name,fill=kName,y=Freq)) +
    geom_bar(stat='identity') +  
    xlab('') +
    ylab('Per cluster (%)') + 
    coord_flip(ylim = c(0,100),xlim=c(0.5,nClusters+0.5),expand = FALSE) +
    myTheme+
    scale_fill_brewer(palette=8,type = 'qual')
  
  gCCN <- ggplot(kTableColN, aes(x=name,fill=kName,y=Freq/1000)) + 
    geom_bar(stat='identity') +  
    xlab('') +
    coord_flip(xlim=c(0.5,nClusters+0.5),expand = FALSE) + 
    myTheme + 
    labs(y=expression(Number~per~cluster~(x10^"3"))) + 
    scale_fill_brewer(palette=8,type = 'qual')
  
  return(list(gFig1=gCKB,
              gFig2=gCC,
              gFig3=gCCN,
              dA=dAll))
}

## Gives count, mean, standard deviation, standard error of the mean, and confidence interval (default 95%).
##   data: a data frame.
##   measurevar: the name of a column that contains the variable to be summariezed
##   groupvars: a vector containing names of columns that contain grouping variables
##   na.rm: a boolean that indicates whether to ignore NA's
##   conf.interval: the percent range of the confidence interval (default is 95%)
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
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
                   c(N      = length2(xx[[col]], na.rm=na.rm),
                     mean   = mean   (xx[[col]], na.rm=na.rm),
                     median = median (xx[[col]], na.rm=na.rm),
                     sd     = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  
  # Rename the "mean" column    
  datac <- rename(datac, c("mean" = measurevar))
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}

#source('scripts/sortingPaper_getH3K4me3_and_H3K9ac_table.R')
#source('~/scripts/R/sortingPaper_getH3K4me3table.R')