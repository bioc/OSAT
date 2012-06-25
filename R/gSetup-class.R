
## Auther: Li Yan
## 

## helper function
## this is one constructor of above the critical part deal with block variable
## assume block ID in sample is 'sFactor', and in wells 'cFactor' 
assign.sample.to.container <- function(smpl, wells, debug=FALSE){
  factorList <- c(colnames(smpl), colnames(wells))
  if (length(factorList) != length(unique(factorList)))
    stop("Common column names exists in sample and container! Cannot go on")
  
  avgCount <- average.2D(smpl, wells, s='sFactor', b='cFactor')
  mCount <- floor(avgCount)     ## interger part
  rCount <- avgCount - mCount   ## residual part
  
  ## step 1, assign mCount
            
  selectedWells <- block.strata.sampling(wells, mCount, 'cFactor', 'sFactor')
  w1 <- selectedWells$selected
  w1 <- w1[order(w1$sFactor, w1$cFactor),]
if(debug) { print("w1");print(head(w1))}

  selectedSamples <- block.strata.sampling(smpl, t(mCount), 'sFactor', 'cFactor')
  s1 <- selectedSamples$selected
  colnames(s1)[3] <- "ID_unit_sample"
  
  ## link selected samples to wells
  part1 <- cbind(smpl[s1[,3],], wells[w1[,3],] )
  if(debug) { print("part1 nrow=");print(nrow(part1))}
  if(debug) print(with(part1, table(cFactor, sFactor))) ## same as mCount
  if(debug) {print(mCount);print(head(part1))}
  
  ## step 2, assign the rest samples
  x <- selectedSamples$residual
  w <- selectedWells$residual
  r <- fraction.strata(x=x, s='sFactor',
                       w=w, b='cFactor', p=rCount)
  if(debug){print("fraction r"); print(head(r))}
  part2 <- cbind.data.frame(x[r[, "sFactor"],], w[r[,"cFactor"],])
  if(debug){ print("part2 nrow=");print(nrow(part2))}
  final <- rbind.data.frame(cbind(part1, stage=1), cbind(part2, stage=2))
}
# tmp <- assign.sample.to.container(smpl, wells, debug=FALSE)
# multi.barplot(tmp, grpVar='plates', varList=c("SampleType", "Race", "AgeGrp"), main="testing assignment")




# TODO: add validation function
get.gSample <- function(x){
  if(!is(x, "gExperimentSetup")) stop("Not a gExperimentSetup class object.\n")
  return(x@data$sample)
}
get.gContainer <- function(x){
  if(!is(x, "gExperimentSetup")) stop("Not a gExperimentSetup class object.\n")
  return(x@data$container)
}


## helper, create expSetup slot
## x: a gExperimentSetup object with @data$sample and @data@container filled
## link: a data frame that indicate the links, have "OrigRow" and "wellID", maybe more
set.expSetup <- function(x, link){
  if(!is(x, "gExperimentSetup")) stop("Not a gExperimentSetup class object.\n")
  # get all columns in sample data frame
  raw <- get.gSample(x)@rawData
  v1 <- colnames(raw)
  rslt <- cbind(raw[link$OrigRow,],
                  link[, !(colnames(link) %in% v1)] )

  # get container columns
  cn <- getLayout(get.gContainer(x))
  v2 <- colnames(cn)
  rslt <- merge(rslt, cn)
  v3 <- colnames(link)[!colnames(link) %in% c(v1, v2)]
  rslt <- rslt[order(rslt$OrigRow), c(v1, v3, v2)]
  rownames(rslt) <- rslt$OrigRow
  return(rslt)
}


setMethod("initialize", "gExperimentSetup",
          function(.Object, sample, container,...){
            if(missing(sample) | missing(container))
              stop("Both sample and container information are needed.\n")
            if(!is(sample,"gSample"))
              stop("Not a gSample object\n")
            if( !is(container, "gContainer"))
              stop("Not a gContainer object\n")
            if(nrow(samples(sample)) > nrow(getLayout(container)) )
              stop("More samples than available wells")
            
            avgCount <- average.2D(sid(sample), cid(container),
                                   s='sFactor', b='cFactor')
            avgCount2 <- average.2D(sid(sample), cid(container),
                                   s='oFactor', b='cFactor')
            link <- assign.sample.to.container(sid(sample), cid(container))
            # per plate block strata counts
            
            blockTable <- table(link$cFactor, link$sFactor)
            optimalTable <- table(link$cFactor, link$oFactor)
            
            .Object@data$sample <- sample
            .Object@data$container <- container
            .Object@data$link <- link
            .Object@expSetup <- set.expSetup(.Object, link)
            .Object@summaryInfo$average.block <- avgCount
            .Object@summaryInfo$average.optimal <- avgCount2
            .Object@summaryInfo$blockTable <- blockTable
            .Object@summaryInfo$optimalTable <- optimalTable
            .Object@metadata$optimalFunction <- "No optimaization."
            .Object
          }
          )
## test
## g1 <-  new("gExperimentSetup", gs, container)

## more friendly constructor
create.experiment.setup <- function(sample, container, ...){
  ## TODO: chech validity
  ##
  object <- new("gExperimentSetup", sample, container)
}

setMethod("samples", "gExperimentSetup", function(x)  return(samples(x@data$sample)))
## samples(gSetup)[1:5,]

expSetup <- function(x){
  if(!is(x, "gExperimentSetup")) stop("Not a gExperimentSetup class object.\n")
  return(x@expSetup)
}

## genetic defined at gContainer
setMethod("get.gAssembly", "gExperimentSetup", function(x) {
  return(x@data$container@data$assembly)}
          )

setMethod("metadata", "gExperimentSetup", function(x)  return(x@metadata)
          )
setReplaceMethod("metadata", "gExperimentSetup",
                 function(x, value) {
                     if (!is.list(value))
                         stop("replacement 'metadata' value must be a list")
                     if (!length(value))
                         names(value) <- NULL # instead of character()
                     x@metadata <- value
                     x
                 })

expLink <- function(x){
  if(!is(x, "gExperimentSetup")) stop("Not a gExperimentSetup class object.\n")
  return(x@data$link)
}

## NOTE: summary already generic, and signature with object, not x
setMethod("summary", "gExperimentSetup", function(object) return(object@summaryInfo ))

# public
get.experiment.setup <- function(x){
  r <- expSetup(x)
  l <- colnames(r) %in% c("OrigRow", "sFactor", "oFactor", "cFactor","stage",
                         "chipID","rowID","wellID")
  return(r[, !l])
}



setMethod("show", "gExperimentSetup", 
function(object) {
        t2 <- summary(object)
	cat("An object of class \"",class(object),"\"\n",sep="")
        #cat("Number of samples by blocking variable\n")
        #print(get.block.strata.table(get.gSample(object)))

        #data <- expSetup(object)
        cat("\nNumber of samples in each plate by blocking strata\n")
        #print(table(data$cFactor, data$sFactor))
        print(t2$blockTable)
        
        #TODO: check if optimal variable exist
        #cat("\nNumber of samples by optimization variable\n")
        #print(get.optimal.strata.table(get.gSample(object)))
        
        cat("\nNumber of samples in each plate by optimization strata\n")
        #print(table(data$cFactor, data$oFactor))
        print(t2$optimalTable)

        x <- get.experiment.setup(object)
	if(length(x) > 0) {
          cat("\nThe experiment setup is", "\n",sep="")
	  cat("@expSetup","\n",sep="")
          if (nrow(x) > 6){
	    print(head(x, 3))
	    cat("\n ... \n")
            print(tail(x, 3))
          }else{
            print(x)
          }
        }
	
})


## map a gExperimentSetup to MSA plate
## rely on the order of wells on the assembly (wellID)
setMethod("map.to.MSA", signature(x="gExperimentSetup", y="MSAroboticPlate"),
 function( x,y) {
   l <- getLayout(get.gAssembly(x))
   filled <- expSetup(x)
   filled$unusedWell <-FALSE
     ## noticed that this is sort by sample order

#   print("nrow in assembly"); print(nrow(l))
#   print("used in setup"); print(nrow(filled))
   out <- merge(l, filled, all.x=TRUE, sort=FALSE)
   out$unusedWell <- ifelse(is.na(out$unusedWell), TRUE, FALSE)
   out <- out[order(out$wellID),]
   idx <- colnames(out) %in% c("chipRows", "chipColumns", "chips", "rows", "columns", "wells",  "chipID", "rowID", "OrigRow","sFactor","oFactor", "stage", "cFactor")
   out <- out[,!idx]
   callGeneric(out, y)
})

# tmp3 <- map.to.MSA(gSetup, MSA4.plate)
# with(tmp3, table(MSA_plate, plates))
# with(tmp3, table(MSA_plate, SampleType))






## assume oFactor exist
## assume cFactor and sFactor
optimal.block <- function(x, nSim=100){
  if(!is(x, "gExperimentSetup")) stop("Not a gExperimentSetup class object.\n")
  #expected counts by group
  oCount <- summary(x)$average.optimal
  optValue <- rep(NA, nSim)  ## hold optimal criteria value
  for (i in 1:nSim) {
    r <- assign.sample.to.container(sid(get.gSample(x)),
                                    cid(get.gContainer(x)))
    
    oCountObs <- as.matrix(table(r$cFactor, r$oFactor))
    
    optValue[i] <-  sum( (oCountObs - oCount)^2)
    if (i==1){
      bestResult <- r
      bestDiff <-  optValue[i]
    }else{
      if (optValue[i] < bestDiff){
        bestResult <- r
        bestDiff <- optValue[i]
      }
    }
  }

  x@data$link <- bestResult
  x@expSetup <- set.expSetup(x, bestResult)
  x@metadata$optimalFunction <- "optimal.block"
  x@metadata$optValue <- optValue
    return(x)
}

# g2 <- optimal.block(g1, nSim=100)









## main function
## sid:    Sample list. A data frame sorted by strata s.
##           each row represent a sample. No NA is allowed in strata variables
## cid:    A data frame represent the experiment wells, sorted by bin b.
##           each row represent a well.
## s:      Strata, a vector of variable names in sid.
##           must be factor or other categorial variable.
## optimalVar:    A vector of variable names in sid,
##         the deviates from frequency are minimized by block and strata.
## b:      Variable indicate block (bins) in the container cid.
## nSim:   number of simulation used in optimization.
##
## return a data frame with variables
##     s optimalVar sid_UNIT b cid_UNIT   ??  iSimulation ?  Prob Stratum
## sid_UNIT:    row number in sid data frame
## cid_UNIT:    row number in cid data frame
##


## new version
optimal.shuffle <- function(x,  nSim=100, k=2){
  if(!is(x, "gExperimentSetup")) stop("Not a gExperimentSetup class object.\n")
  oCount <- summary(x)$average.optimal
  con <- getLayout(get.gContainer(x))

  link <- expLink(x)
    NAMES <- colnames(link)
    v1 <- NAMES[!NAMES %in% colnames(con)] # sample related columns
    v2 <- NAMES[NAMES %in% colnames(con)]  # container related columns
  ## observed counts  
  oCountObs <- as.matrix(table(link$cFactor, link$oFactor))
  
  optValue <- rep(NA, nSim)  ## hold optimal criteria value
  bestResult <- link
  bestDiff <- sum( (oCountObs - oCount) ^ 2 )
  optValue[1] <- bestDiff
  for (i in 2:nSim) {
        ## shuffle
    s <- sample(link[,'wellID'],k)       # sample k wells
      idx <- which(link$wellID %in% s)
    dd <- link[idx,]                     # get the k records
    d1 <- dd[,v1]                        # sample columns
    d2 <- dd[,v2]                        # container columns
      idx1 <- (2:(k+1) ) %%k
      idx1[which(idx1==0)] <- k          
    newdd <- cbind(d1, d2[idx1,])        # shuffle container (wells)
    #newdd <- newdd[, colnames(link)]    # keep same order of columns
    newlink <- rbind(link[-idx,], newdd)
    
    oCountObs <- as.matrix(table(newlink$cFactor, newlink$oFactor))
    optValue[i] <- sum( (oCountObs - oCount) ^ 2 )
    if (optValue[i] < bestDiff){
        bestResult <- newlink
        bestDiff <- optValue[i]
        link <- newlink
      }
  }
  
  x@data$link <- bestResult
  x@expSetup <- set.expSetup(x, bestResult)
  x@metadata$optimalFunction <- "optimal.shuffle"
  x@metadata$optValue <- optValue
    return(x)  
                 
}

## g3 <- optimal.shuffle(g1, nSim=10000, k=2)

## TODO: passing function as parameters
## FUN <- match.fun(function_name); FUN(data)
## or f=function(x, f), eval(as.call(list(as.name(f), x)))


create.optimized.setup <- function(fun="default", sample, container, ...){
  #TODO: check validity of param
  .Object <- new("gExperimentSetup", sample, container)
  param <- list(...)
  param$x <- .Object
  
  ## dispatch to different object function
  if (missing(fun) | fun=="default"){
    warning("Using default optimization method: optimal.shuffle")
    .Object <- do.call("optimal.shuffle", param)
  } else{
    .Object <- do.call(fun, param)
  }
}

## g4 <- create.optimized.setup(sample=gs, container=container, nSim=10000, k=2)
## g5 <- create.optimized.setup("optimal.shuffle", sample=gs, container=container, nSim=10000, k=2)
## g6 <- create.optimized.setup("optimal.block", sample=gs, container=container, nSim=100)


## ## TODO: turn to a class
## f1 <- function(x, n){rep(x, n)}
## f1('d', 3)
## list1 <- list(FUN='f1', n=4, x='d') ## other parameters following
## ff1 <- function(o) {
##   FUN <- match.fun(o$FUN)
##   FUN(o$x, o$n)
## }
## ff1(list1)

## list2 <- list(n=5, x='d')
## ff2 <- function(fun, o){
##   do.call(fun, o)
## }
## ff2('f1', list2)



## make a toal random design

create.random.setup <- function(sample, container, ...){
  x <- new("gExperimentSetup", sample, container)
  random.assign(x, ...)
}

## x is a gSetup object
## r is a (random) vector length nrow(container)
random.assign <- function(x, r){
  #print("To select wells randomly and link to sample")
  if(!is(x, "gExperimentSetup")) stop("Not a gExperimentSetup class object.\n")
  w <- getLayout(get.gContainer(x))
  if (missing(r))  r <- runif(nrow(w))
  if (nrow(w) != length(r))
    stop ("Length of random vector not equal to available container wells")

  w <- w[order(r),]
  
  s <- sid(get.gSample(x))
  
  # print(nrow(w));print(head(w));
  a0 <- cbind(s, w[1:nrow(s),], stage=0)  # linked
  # print(head(a0))
  link <- a0[,colnames(expLink(x))]  # replace link in setup

  oCount <- summary(x)$average.optimal
  oCountObs <- as.matrix(table(link$cFactor, link$oFactor))
  optValue <- sum( (oCountObs - oCount) ^ 2 )
 
  x@data$link <- link
  x@expSetup <- set.expSetup(x, link)
  x@metadata$optimalFunction <- "random"
  x@metadata$optValue <- optValue
    return(x)                 
}



# tmp <- random.assign(gSetup)
# QC(tmp)

# tmp1 <- create.random.setup(gs, gc)

## r <- sample(nrow(getLayout(gc)), nrow(getLayout(gc)), replace=FALSE)
## tmp2 <- create.random.setup(gs, gc, r)
## QC(tmp2)
