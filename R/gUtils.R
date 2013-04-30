
## randomlly assign samples to container
# samples: a data frame representing list of samples
# wells:   a data frame representing container
#          nrow(samples) <= nrow(wells)
random.setup <- function(samples, wells){
    c0 <- cbind(wells, r=runif(nrow(wells)))
    c0 <- c0[order(c0$r),]
    a0 <- cbind(samples, c0[1:nrow(samples),])
    return(a0)
}

## multiple barplot in one page
# x:       a data frame, count of rows by grpZVar will be the barplot y-axis
# grpVar:  a variable name in x, barplot will 
# varList: a vector of variable names in x, one plot per elements 
# main:    the main text on of the barplot
multi.barplot <- function(x, grpVar="plates", varList, main=NULL, ...){
    nPlot <- length(varList)
    nPlotCol <- ceiling(sqrt(nPlot))
    nPlotRow <- ceiling(nPlot / nPlotCol)

    nPlates <- length(unique(x[, grpVar]))

    layout(matrix(1:(nPlotCol*nPlotRow), nPlotRow, nPlotCol, byrow = TRUE))
    for (i in varList){
      cnt <- table(x[,c(grpVar, i)])
      nLevels <- length(dimnames(cnt)[[2]])
      colScale <- rep(seq(0.1,0.9,length=nPlates), nLevels)
      ## TODO: pass part of ... to here
      barplot(cnt,  xlab=i, col=gray(colScale),  beside=TRUE, ...)
    }
    if (is.null(main)) main="Sample distribution"
    mtext(main, side=3, outer=TRUE, line=-3) 
}


# http://svn.r-project.org/R/trunk/src/library/stats/R/htest.R

## do Chisq-test compare counts in grpVar, on each variables in varList 
multi.chisq.test <- function(x, grpVar="plates", varList, main=NULL){
    DNAME <- deparse(substitute(x))
    rslt <- NULL
    for (i in varList){
      cnt <- table(x[,c(grpVar, i)])
      t <- chisq.test(cnt)
      t1 <- cbind.data.frame(i,
                       statistic=t$statistic,
                       parameter=t$parameter,
                       p.value=t$p.value)
      rslt <- rbind(rslt, t1)
      #print(chisq.test(cnt))
    }
    colnames(rslt) <- c("Var", names(t$statistic),  names(t$parameter),  "p.value")
    rownames(rslt) <- 1:nrow(rslt)
    return(list(method=t$method, data.name=DNAME, stat=rslt))
}


## average count of each sample strata per container bin
  ## sid:    sample, strata by one variable s
  ## cid:    container, blocked by one variable b
  ## return matrix:
  ##   each row represent average count per container block (ex: chips)
  ##   each column represent a level in sample strata (ex: Race)
  ##   colnames and rownames are strata and block factors, respectively
average.2D <- function(sid, cid, s, b){
  if ( length(s)!=1 | length(b) !=1 )
    stop("Cannot handle more than one variables in each data frame")
  optCount <- as.data.frame(table(sid[,s], useNA = "no"))## count per strata
  binList <- as.data.frame(table(cid[,b], useNA = "no")) ## count per block
  ## calculate average counts per group per chip
  f <- function(x) x*binList$Freq/sum(binList$Freq)
  optCount1 <- sapply(optCount$Freq,f)
  colnames(optCount1) <- optCount[,1]
  rownames(optCount1) <- binList[, 1]
  return(optCount1)
}

# remove all zero rows and columns
remove.zeros <- function(x){ #removel rows and columns with all 0
    idx1 <- !(apply(x, 1, function(x) all(x==0)) )
    idx2 <- !(apply(x, 2, function(x) all(x==0)) )
    x <- x[idx1, idx2]  
}





## new function
  ## x: a data frame to be sampled, sorted and stratified by varaible s
  ## m: an integer matrix represent multiple non-replacement sampling from x
  ##    each row represent the sample sizes for each strata i, in which
  ##    each column represent one sampling from x,
  ##    clearly we need sum(m) <= nrow(x) and rowSums(m)[s] <= n(x) in each s
  ##    its row name should corresponding to factors in s
  ##    its colnames should corresponding to 
  ## s: a variable name that indicate the strata in x, corresponding the rows in m
  ## b: a variable name that corresponding to columns in m
  ## return a list
  ##    selected: samples selected from x
  ##    residual: the residual in x after selection
block.strata.sampling <- function(x, m, s, b){
  dataRow <- 1:nrow(x)
  selected <- NULL
  for ( i in 1:nrow(m)){                   # strata loop
    y <- NULL
    iIDX <- rownames(m)[i]                 # strata name
    d <- dataRow[x[,s] == iIDX]            # data in the strata
    ni=sum(m[i,])
    if(ni>0){
      y <- sample(d, ni, replace = FALSE)  # select from strata
        #  split to block (corresponding to columns in m)
      yy <- cbind.data.frame(rep(iIDX, ),
                             rep(colnames(m), m[i,]), y)
      selected <- rbind(selected, yy)
    }
    
  }
  colnames(selected) <- c(s, b, "ID_unit")
  return(list(selected=selected, residual=x[-selected[,3],]))
}

## set.seed(123)
## tmp <- block.strata.sampling(wells, m, "cFactor", "sFactor")

## set.seed(123)
## tmp1 <- block.strata.sampling(wells, m, "cFactor", "sFactor")
## tmp2 <- block.strata.sampling(smpl, t(m), "sFactor", "cFactor")

## identical(tmp$residual, tmp1$residual)





## distribute samples into blocked (clustered) bins by strata with unequal probability
##
  ## x: sample data frame, strata by variable s='combinedIdx'
  ## s: a variable represent strata in sample (factor)
  ## w: container data frame, each row is a bin,
  ##    bins are blocked by variable b= chipID (blockIdxName)
  ## b: a variable represent blocks in container
  ## p: a matrix of inclusion probabilities or auxiliary information
  ##      used to compute them. It represent probability of sample
  ##      assigned into a block of  bin. 
  ##    each column represent one strata from x,
  ##      ncol(p) == length(unique(x$s))
  ##    each row represent the probobility the block of bin is chosen
  ##      nrow(p)==length(unique(w$b)
  ##    probability do not need to be normalized to one
  ##    its rowname should corresponding to strata in s
  ##    its colnames should corresponding to levels in w
fraction.strata <- function(x, s='sFactor', w, b='cFactor', p){
  rslt <- NULL
  colIdx <- colnames(p)
  xIdx <- 1:nrow(x)
  wIdx <- 1:nrow(w)     # all wells
  wSelected <- NULL
  for (i in colIdx){        # each col represent a sample strata
    if (all(p[,i]==0)) next # no sample need to assign
    prob <- p[,i, drop=FALSE] 
    w1 <- wIdx[w[,b] %in% rownames(prob)[prob==0]] # selection prob ==0
    w2 <- wIdx[as.vector(wSelected)]   # previously selected

    if (length(c(w1, w2)) != 0){
      n <- wIdx[-c(w1, w2)]  # available wells
    }else{
      n <- wIdx
    }

    nw <- as.data.frame(table(w[n,b]))  # count of wells per bin
    ww <- prob/nw[, 'Freq']
    # expand to each well. sapply can handle NaN, Inf, and 0
    pik <- sapply(1:nrow(ww), function(j) ww[rep(j, each = nw$Freq[j]), ], simplify = FALSE)
    pik <- unlist(pik)
    # must prevent n to be treated as number
    tIdx <- as.numeric(as.vector(sample(as.character(n), length(xIdx[x[,s]==i]), replace = FALSE, prob = pik)))
    rslt <- rbind(rslt, cbind(xIdx[x[,s]==i], tIdx) )
    wSelected <- c(wSelected, tIdx)
  }
  colnames(rslt) <- c(s, b)
  return(rslt)
}



barplot.gExperimentSetup <- function(object, main=NULL, ...){
    data <-  get.experiment.setup(object)
    optValues <- metadata(object)$optValue
    
    sVarList <- get.gSample(object)@strata
    oVarList <- get.gSample(object)@optimal
    varList <- c(sVarList, oVarList[!oVarList %in% sVarList])

    if (!is.null( optValues)){
      nPlot <- length(varList) +1  # add obj value plot
    }else{
      nPlot <- length(varList)
    }
    
    nPlotCol <- ceiling(sqrt(nPlot))
    nPlotRow <- ceiling(nPlot / nPlotCol)

    blockLevel <- get.gContainer(object)@batch
    nPlates <- get.gContainer(object)@n

    
    #match.arg(..., "layout")
    layout(matrix(1:(nPlotCol*nPlotRow), nPlotRow, nPlotCol, byrow = TRUE))
    for (i in varList){
      cnt <- table(data[,c(blockLevel, i)])
      nLevels <- length(dimnames(cnt)[[2]])
      colScale <- rep(seq(0.1,0.9,length=nPlates), nLevels)
      #print(colScale)
      ## TODO: pass part of ... to here
      barplot(cnt,  xlab=i, col=gray(colScale),  beside=TRUE)
    }

    # plot the opttmization objective function values
    
    if (!is.null( optValues)){
    plot(optValues, ylab="Value of objective function")
    x <- which( optValues  == min(optValues))[1]
    y <- min(optValues)
    points(x, y, pch=23, col="red",  cex=2.5)
    points(1, metadata(object)$optValue[1], pch=23,col="blue", cex=2.5)
    }
    
    if (is.null(main)) main="Sample distribution"
    mtext(main, side=3, outer=TRUE, line=-3) 
}





setMethod("plot",c("gExperimentSetup","missing"),
          function(x,y,main=NULL, ...){
          barplot.gExperimentSetup(x, main, ...)
})


QC <- function(object, main=NULL, ...){
    if (class(object) != "gExperimentSetup"){
      stop("Only for object in gExperimentSetup class")
    }
    data <-  get.experiment.setup(object)

    sVarList <- get.gSample(object)@strata
    oVarList <- get.gSample(object)@optimal
    varList <- c(sVarList, oVarList[!oVarList %in% sVarList])
    
    blockLevel <- get.gContainer(object)@batch
    plot(object, grpVar=blockLevel, varList=varList, main=main)

    t <- invisible(multi.chisq.test(data, grpVar=blockLevel, varList=varList))
    message("\nTest independence between \"", blockLevel, "\" and sample variables\n", sep="")
    message("\n", t$method, "\n")
    print(t$stat)
    message("\n")
}

# save(BeadChip96ToMSA4MAP, GenotypingChip , IlluminaBeadChip, IlluminaBeadChip24Plate, IlluminaBeadChip48Plate, IlluminaBeadChip96Plate, MSA4.plate, file="predefined.Rdata")
