
# Author: Li Yan


## helper function
## combine multiple factor columns to one
## TODO: see strata in "sampling" for different way (unique() function)
combineFactors <- function(x, s, outNames=c('origRow', 'newFactor')){
      cnt <- as.data.frame(table(x[,s, drop=FALSE]), useNA = "ifany")
      if (length(s)==1) colnames(cnt)[1] <- s
      cnt$newFactorLevel <- as.factor(1:nrow(cnt))

      tmp <- cbind(x[,s, drop=FALSE], OriginalRowIdx = 1:nrow(x))
      a <- merge(tmp, cnt, by=s) 
      a <- a[order(a$OriginalRowIdx),c(s, 'OriginalRowIdx', 'newFactorLevel')]
      colnames(a) <- c(s, outNames)
      colnames(cnt)[ncol(cnt)] <- outNames[2]
      return(list(data=a, cnt=cnt))
  }

#### method for gSample class

## class constructor
## change a data frame to class
## return: an gSample object
setMethod("initialize", "gSample",
          function(.Object, x, optimal, strata, ...){
            vList <- !(colnames(x) %in% unique(c(strata, optimal)) )
            #if (!is.null(vList))
              annotation <- x[,  vList]
            smpl <- combineFactors(x, strata, outNames=c('OrigRow', 'sFactor'))
            opt <- combineFactors(x, optimal, outNames=c('OrigRow', 'oFactor'))
            sid=merge(smpl$data, opt$data,
              by=c('OrigRow', intersect(strata, optimal)), all=TRUE)
            sid <- sid[order(sid[,'OrigRow']),]
            # only keep minimum info
            rownames(sid) <- sid$OrigRow
            sid <- sid[,c("OrigRow", "sFactor", "oFactor")]
            
            .Object@rawData <- x
            .Object@optimal <- optimal
            .Object@strata <- strata
            .Object@data$optimalTable <- opt$cnt
            .Object@data$strataTable <- smpl$cnt
            .Object@data$annotation <- annotation
            .Object@data$sid <- sid
            .Object
          }
          )

## user friendly class constructor
gSample <- function(x, optimal, strata){
  if (missing(strata))  strata=optimal[1]
  object <- new("gSample", x=x, optimal=optimal, strata=strata)
  return(object)
}

## public
## user friendly way of construct a gSample object
setup.sample <- function(x, optimal, strata ){
  ## if ( class(x)[1] == "ExpressionSet") {xx <-  pData(x)
  ##                                     }else{
  ##                                       xx <- x}
  ## gSample(xx, optimal,strata )
  gSample(x, optimal,strata )
}

## setMethod("setup.sample", "data.frame", function(x, optimal, strata ){
##   gSample(x, optimal,strata )
## })
## setMethod("setup.sample", "ExpressionSet", function(x, optimal, strata ){
##   x <- pData(x)
##   gSample(x, optimal,strata )
## })

## DO not create setter to avoid inconsistency

## some getters


setMethod("samples", "gSample", function(x)  return(x@rawData))
setMethod("sid", "gSample", function(x)  return(x@data$sid))

setMethod("summary", "gSample", function(object) 
  return (list(strataTable=object@data$strataTable,
               optimalTable=object@data$optimalTable)))

setMethod("show","gSample",
#  Print and show method 
function(object) {
	cat("An object of class \"",class(object),"\"\n",sep="")
	x <- slot(object,"rawData")
        d <- summary(object)
	if(length(x) > 0) {
          cat("The raw data are", "\n",sep="")
	  if (nrow(x) > 6){
	    print(head(x))
	    cat("\n ... \n")
            print(tail(x))
          }else{
            print(x)
          }
          cat("\n Blocking strata in the data:\n\n",sep="")
          print(d$strataTable)
          cat("\n Optimization strata in the data\n\n",sep="")
          print(d$optimalTable)
        }
	
})

## test
## gs <- setup.sample(pheno,  optimal=c("SampleType", "Race", "AgeGrp"), strata=c("SampleType") )
## gsTmp <- setup.sample(pheno,  optimal=c("SampleType", "Race", "AgeGrp"))
## identical(gs, gsTmp)




## TODO, directly from txt file
## read.sample <- function(file, ...){
## }


