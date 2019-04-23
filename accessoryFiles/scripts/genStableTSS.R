args <- commandArgs(trailingOnly = TRUE)
nFold <- as.numeric(args[1])
fIn   <- args[2]
fOut  <- args[3]

#r <- read.table('allTSS.txt')
r <- read.table(fIn)

## Set parameters
nLow  <- 0.15 ## exclude bottom 15%
nHi   <- 0.01 ## exclude top 1%
#nFold <- 1.40 ## require <X-fold in all samples

## make Generic names
sNames <- paste0("s",1:(dim(r)[2]-3))
names(r) <- c("cs","from","to",sNames)

# Remove below BG signals
hTSS <- r
for (n in sNames){
  hTSS <- hTSS[hTSS[,n] > 0,] 
  print (dim(hTSS))
}

## Exclude out-of-boundry TSSs
hTSS$ok  <- TRUE

for (n in sNames){
  lo  <- quantile(hTSS[[n]],nLow,na.rm=TRUE)
  hi  <- quantile(hTSS[[n]],1-nHi,na.rm=TRUE)
  hTSS$ok[ (hTSS[[n]] >= hi)  | (hTSS[[n]] <= lo)]   <- FALSE
}

## Make "final" set of clean TSSs
okTSS <- hTSS[hTSS$ok, ]
  
## ALL v ALL comparison
okTSS$allOK <- TRUE
  
for (n1 in sNames){
  for (n2 in sNames){
    ## Sum of signals
    okTSS[,paste0('tot_',n1,'_v_',n2)] <- (okTSS[,n1]+okTSS[,n2])/2
    
    ## Log2 difference - median(log2 diff)
    okTSS[,paste0('l2_',n1,'_v_',n2)]  <- log2(okTSS[,n1]/okTSS[,n2])  - log2(median(okTSS[,n1]/okTSS[,n2]))
    
    ## replace NaNs with 0s
    okTSS[is.nan(okTSS[,paste0('l2_',n1,'_v_',n2)]),paste0('l2_',n1,'_v_',n2)] <- 0
    
    ## Mark TSS with log2(delta) > log2(nFold)
    okTSS$allOK[abs(okTSS[,paste0('l2_',n1,'_v_',n2)]) >= log2(nFold)]  <- FALSE
  }
}

## Get stable TSSs
stableRegions <- okTSS[okTSS$allOK,1:(dim(r)[2])]

## write to output file
write.table(x = stableRegions[,1:3], file = fOut, 
            quote = FALSE, sep = "\t", row.names = FALSE,
            col.names = FALSE)