## define classes representing containers (chip, plate etc)
## used in genomic experiments, and classes representing the assembly
## of all plate/chips/wells
## Auther: Li Yan
## 

#### methods for gArray class
## getter and setter of metadata, use with care
setMethod("metadata", "gArray",
          function(x) {
              if (is.null(x@metadata) || is.character(x@metadata))
                  list(metadata = x@metadata)
              else
                  x@metadata
          })
setReplaceMethod("metadata", "gArray",
                 function(x, value) {
                     if (!is.list(value))
                         stop("replacement 'metadata' value must be a list")
                     if (!length(value))
                         names(value) <- NULL # instead of character()
                     x@metadata <- value
                     x
                 })



setMethod("nrow", "gArray", function(x) x@nRows)
setMethod("ncol", "gArray", function(x) x@nColumns)
setMethod("dim", "gArray",  function(x) cbind(nrow(x), ncol(x)))
setMethod("getLayout", "gArray", function(x, ...) x@layout)



#### methods for gSlide class
setMethod("initialize", "gSlide",
          function(.Object, nRows, nColumns, byrow, comment, ...){
            if(!missing(nRows))     .Object@nRows <- as.integer(nRows)
            if(!missing(nColumns))  .Object@nColumns <- as.integer(nColumns)
            if(!missing(byrow))     .Object@byrow <- byrow
            if(!missing(comment))   .Object@metadata$comment <- comment
            
            if (.Object@byrow){
              .Object@layout <-
                data.frame(matrix(c(rep(1:.Object@nRows, each=.Object@nColumns),
                                    rep(1:.Object@nColumns,.Object@nRows)), ncol=2))
            } else{
              .Object@layout <-
                data.frame(matrix(c(rep(1:.Object@nRows, .Object@nColumns),
                                    rep(1:.Object@nColumns,each=.Object@nRows)), ncol=2))
            }
            .Object@layout <- cbind(.Object@layout, 1:nrow(.Object@layout))
            rownames(.Object@layout) <- 1:nrow(.Object@layout)
            colnames(.Object@layout) <- c("rows", "columns", "wells")
            .Object
          }
          )

#  Print and show method 
setMethod("show","gSlide",
function(object) {
	cat("An object of class \"",class(object),"\"\n",sep="")
        cat(metadata(object)$comment, "\n",sep="")
        
	x <- slot(object,"layout")
	if(length(x) > 0) {
          cat("It has", nrow(object), "rows and", ncol(object),"columns.\n", sep=" ")
          cat("The layout is \n@layout\n",sep="")
	  print(x)
	  cat("\n")
        }
	
})

## TODO: constructor to keep the name of object created


#### methods for Beadchip class
## define know chip configurations (as object)
IlluminaBeadChip <- new("BeadChip",  nRows=6, nColumns=2,
                        comment="Illumina Bead Chip have 6 rows and 2 columns.")

GenotypingChip <- new("BeadChip",  nRows=12L, nColumns=1L,
                        comment="Genotyping Chip have 12 rows and 1 columns.")


#### methods for gPlate class
setMethod("initialize", "gPlate",
          function(.Object, chip, nRows, nColumns, comment, ...){
            if(!missing(chip)){    .Object@chip <- chip
                                  }else{
                                    .Object@chip <- IlluminaBeadChip
                                  }
            if(!missing(nRows)){     .Object@nRows <- nRows
                                   }else{
                                     .Object@nRows=2L
                                   }
            if(!missing(nColumns)){  .Object@nColumns <- nColumns
                                     .Object@nColumns <- 4L
                                     }
            if(!missing(comment))   .Object@metadata$comment <- comment

            .Object@metadata$chipName <- deparse(substitute(chip))

            aChip <- .Object@chip
            m <- aChip@layout
            nr <- .Object@nRows
            nc <- .Object@nColumns
            rList <- rep(1:nr, each=nc*nrow(m))  ## chipRow
            cList <- rep(rep(1:nc, each=nrow(m)), nr)
            .Object@layout <- m[rep(seq_len(nrow(m)),nr*nc),]
            .Object@layout <- cbind(rList,cList, rep(1:(nr*nc), each=nrow(m)),
                                    .Object@layout)
            colnames(.Object@layout)[1:3] <- c("chipRows","chipColumns","chips")
            rownames(.Object@layout) <- 1:nrow(.Object@layout)
            .Object
          }
          )

setMethod("show","gPlate",
#  Print and show method 
function(object) {
        nameOfchip <- metadata(object)$chipName   #deparse(substitute(object@chip))
 	cat("\nAn object of class \"",class(object),"\"\n",sep="")
        cat(metadata(object)$comment, "\n",sep="")
        cat("It has", nrow(object), "rows and", ncol(object),"columns of",
            nameOfchip,"chips.\n", sep=" ")
        #TODO: info about chip
	x <- slot(object,"layout")
	if(length(x) > 0) {
          cat("The layout is", "\n",sep="")
	  cat("@layout","\n",sep="")
          if (nrow(x) > 12){
	    print(head(x))
	    cat("\n ... \n")
            print(tail(x))
          }else{
            print(x)
          }
        }
	
})

#### methods for BeadPlate class
## define know plate configurationes
IlluminaBeadChip96Plate <- new("BeadPlate", chip=IlluminaBeadChip, nRows=2L, nColumns=4L, comment="Illumina BeadChip Plate with 8 chips and 96 wells")
IlluminaBeadChip48Plate <- new("BeadPlate", chip=IlluminaBeadChip, nRows=1L, nColumns=4L, comment="Illumina BeadChip Plate with 4 chips and 48 wells")
IlluminaBeadChip24Plate <- new("BeadPlate", chip=IlluminaBeadChip, nRows=1L, nColumns=2L, comment="Illumina BeadChip Plate with 2 chips and 24 wells")


#### methods for assembly class
## TODO: diff way to find the plate name
setMethod("initialize", "gAssembly",
          function(.Object, plate, n, comment){
            if(!missing(plate)){     .Object@plate <- plate
                                     }else{
                                       .Object@plate <- IlluminaBeadChip96Plate
                                       }
            if(!missing(n)){     .Object@n <- as.integer(n)
                               }else{
                                 .Object@n <- 2L
                                 }
            if(!missing(comment))   .Object@metadata$comment <- comment

            .Object@metadata$plateName <- deparse(substitute(plate))
            
            aPlate <- .Object@plate
            m <- getLayout(aPlate)
            n <- .Object@n

            .Object@layout <- m[rep(seq_len(nrow(m)),n),]
            .Object@layout <- cbind(rep(1:n, each=nrow(m)),
                                    .Object@layout)
            colnames(.Object@layout)[1] <- c("plates")
            rownames(.Object@layout) <- 1:nrow(.Object@layout)

            container <- .Object@layout
            nc=length(unique(m[,"chips"]))   # number of chips per plate
            nr=length(unique(m[,"rows"]))    # number of rows on each plate
            chipID <- container[,"chips"]+(container[,"plates"] -1)*nc
            rowID <- container[,3]+ (chipID-1)*nr
            wellID <- 1:nrow(container)
            
            .Object@layout <- cbind(.Object@layout,chipID, rowID, wellID)
            .Object
          }
          )


setMethod("metadata", "gAssembly",
          function(x) {
              if (is.null(x@metadata) || is.character(x@metadata))
                  list(metadata = x@metadata)
              else
                  x@metadata
          })
setReplaceMethod("metadata", "gAssembly",
                 function(x, value) {
                     if (!is.list(value))
                         stop("replacement 'metadata' value must be a list")
                     if (!length(value))
                         names(value) <- NULL # instead of character()
                     x@metadata <- value
                     x
                 })

setMethod("show","gAssembly",
#  Print and show method 
function(object) {
	cat("An object of class \"",class(object),"\"\n",sep="")
        cat(metadata(object)$comment, "\n",sep="")
        cat("It consists of", object@n, metadata(object)$plateName, "plates.\n")
        #TODO: info about plate and chip
	x <- slot(object,"layout")
	if(length(x) > 0) {
          cat("The layout is", "\n",sep="")
	  cat("@layout","\n",sep="")
          if (nrow(x) > 12){
	    print(head(x))
	    cat("\n ... \n")
            print(tail(x))
          }else{
            print(x)
          }
        }
})


setMethod("getLayout", "gAssembly", function(x, ...) x@layout)



#### methods for gContainer class
## TODO: need improvement
# transfer between factor and its value, see ?fact
# as.numeric(levels(f))[f]
## TODO: rowname of the container?
setMethod("initialize", "gContainer",
          function(.Object, plate, n, batch, exclude){
            if(missing(plate) | missing(n) ){
              stop("Need both plate configuration and number.")
            }
            if (batch=="chips" | missing(batch)) {
              blockIdxName <- "chipID"
            } else if (batch=="plates") {
              blockIdxName <- "plates"
            }else{
              stop("Only 'chips' or 'plates' allowed in block level.")
            }
            names(batch) <- blockIdxName   # this is the column name in assmbly
            
            if (missing(exclude)) exclude <- data.frame()

            a <- new("gAssembly", plate, n, comment="new Assembly")
            .Object@data$assembly <- a

            container=getLayout(a)
            container$cFactor <- factor(container[,blockIdxName])

            if ( nrow(exclude)>0){
              ## TODO: use setdiff()
              e <- exclude[,colnames(exclude) %in%
                           c("plates", "chips", "wells"), drop=FALSE]
              e$inExclude <- TRUE
              t <- merge(container, e, sort=FALSE, all=TRUE)
              container <- t[is.na(t$inExclude),!(colnames(t) %in% c("inExclude"))]
            }

             cnt <- as.data.frame(table(container[,"cFactor", drop=FALSE]))
             colnames(cnt)[1] <- "cFactor"
             cnt <- cbind(as.numeric(levels(cnt[,1]))[cnt[,1]], cnt)
             colnames(cnt)[1] <- blockIdxName
                      
             cid <- container[,c(blockIdxName, "wellID", "cFactor")]
             cid <- cid[order(cid$wellID),]

            .Object@plate <- plate
            .Object@n <- as.integer(n)
            .Object@batch <- batch
            .Object@exclude <- exclude
            .Object@data$container <- container
            .Object@data$cid <- cid
            .Object@data$blockTable <- cnt
            .Object
          }
          )
setMethod("getLayout", "gContainer", function(x, ...) x@data$container)
setMethod("metadata", "gContainer",
          function(x) {
              if (is.null(x@metadata) || is.character(x@metadata))
                  list(metadata = x@metadata)
              else
                  x@metadata
          })
setReplaceMethod("metadata", "gContainer",
                 function(x, value) {
                     if (!is.list(value))
                         stop("replacement 'metadata' value must be a list")
                     if (!length(value))
                         names(value) <- NULL # instead of character()
                     x@metadata <- value
                     x
                 })

## get.batch.info <- function(x){
##   if(!is(x, "gContainer")) stop("Not a class gContainer object.\n")
##   return(x@batch)
## }
## use getLayout()
## get.container <- function(x){
##   if(!is(x, "gContainer")) stop("Not a class gContainer object.\n")
##   return(x@data@container)
## }


setMethod("getLayout", "gContainer", function(x)  return(x@data$container))
setMethod("cid", "gContainer", function(x)  return(x@data$cid))
setMethod("get.gAssembly", "gContainer", function(x)  return(x@data$assembly))

## cid(container)[1:5,]


## TODO: more metadata? move plate/n/ etc to metadata and move container out?
setMethod("show", "gContainer", 
function(object) {
	cat("An object of class \"",class(object),"\"\n",sep="")
        cat(metadata(object)$comment, "\n",sep="")
        cat("It consists of", object@n, metadata(object)$plateName, "plates.\n")
        if(nrow(object@exclude)>0){
          cat("There are", nrow(object@exclude), "wells excluded.\n")
        }
        cat("The block level is set at", object@batch, "level.\n")
        print(object@data$blockTable)
 
	x <- object@data$container
	if(length(x) > 0) {
          cat("The container layout is", "\n",sep="")
	  cat("@data$container","\n",sep="")
          if (nrow(x) > 12){
	    print(head(x))
	    cat("\n ... \n")
            print(tail(x))
          }else{
            print(x)
          }
        }
	
})



## TODO: add availibility indicator?

#class() is(c4, "gSlide") is.numeric(c4@layout)
gContainer <- function(plate, n, batch=c("chips", "plates"), exclude=NULL, ...){
  if(missing(plate) | missing(n) ){
    stop("Need both plate configuration and number.")
  }
  if(is.null(exclude)){
    obj <- new("gContainer", plate, n, batch=batch)
  }else{
    obj <- new("gContainer", plate, n, batch=batch, exclude)
  }
  ## a patch ## TODO: move to gAssembly
  metadata(obj)$plateName <- deparse(substitute(plate))
  return(obj)
}

# container <- setup.container(IlluminaBeadChip96Plate, 6, batch='plates')

# a wrapper
setup.container <- function(plate, n, batch="plates", exclude=NULL){
  obj <- gContainer(plate, n, batch, exclude=NULL)
  metadata(obj)$plateName <- deparse(substitute(plate))
  return(obj)
}


## TODO
removeWells <- function(container, wells, ...){
}

assignWells <- function(container, wells, ...){
}





#### methods for MSAroboticPlate class

setMethod("initialize", "MSAroboticPlate",
          function(.Object, nRows, nColumns, byrow, comment, ...){
            if(!missing(nRows))     .Object@nRows <- as.integer(nRows)
            if(!missing(nColumns))  .Object@nColumns <- as.integer(nColumns)
            if(!missing(byrow))     .Object@byrow <- byrow
            if(!missing(comment))   .Object@metadata$comment <- comment

            if (.Object@byrow){
              n1=.Object@nRows
              n2=.Object@nColumns
            }else{
              n1=.Object@nColumns
              n2=.Object@nRows
            }

            .Object@wells <- paste(rep(LETTERS[1:n2], n1),sprintf("%02d",rep(1:n1, each=n2)), sep="")
            .Object@layout <- data.frame(matrix(paste(rep(LETTERS[1:n2], n1),sprintf("%02d",rep(1:n1, each=n2)), sep=""), ncol=n1))
             colnames(.Object@layout) <- 1:n1
             rownames(.Object@layout) <- LETTERS[1:n2]

             .Object
          }
          )

#  Print and show method 
setMethod("show","MSAroboticPlate",
function(object) {
	cat("An object of class \"",class(object),"\"\n",sep="")
        cat(metadata(object)$comment, "\n",sep="")
        
	x <- getLayout(object)
	if(length(x) > 0) {
          cat("It has", nrow(object), "rows and", ncol(object),"columns.\n", sep=" ")
          cat("The layout is \n@layout\n",sep="")
	  print(x)
	  cat("\n")
        }
	
})

# tmp <- new("MSAroboticPlate")

## no need to define generic yet, only one know robotic plate so far
wells <- function(x){
  if(!isClass(x, "MSAroboticPlate")) stop("Not a MSAroboticPlate class object")
  return(x@wells)
}

MSA4.plate <- new("MSAroboticPlate", nRows=8L, nColumns=12L, byrow=FALSE,
           comment="MSA-4 robotic plate")



setMethod("map.to.MSA", signature(x="data.frame", y="MSAroboticPlate"),
 function(x, y) {
   # number of MSA plate needed
   w <- wells(y)
   nw <- length(w)
   n <- ceiling(nrow(x)/nw)
   
   ## # a list of plate positions
   l <- cbind.data.frame(MSA_plate=rep(1:n, each=nw),
                         MSA_ID=rep(w, n))

   nn <- n*nw - nrow(x) # empty rows
   if (nn>0) {
     xx <- as.data.frame(matrix(rep(NA, nn*ncol(x)), ncol=ncol(x) ))
     colnames(xx) <- colnames(x)
     x <- rbind(x, xx)
   }
   
   return(cbind(MSA_plate=rep(1:n, each=nw),
                MSA_ID=rep(w, n),
                x))
   
})

# tmp1 <- map.to.MSA(get.expSetup(gSetup), MSA4.plate)  ## wrong chip order

## map BeadChip Plate to robotic plate
## notied that the 'chips' is the "BC position"
BeadChip96ToMSA4MAP <- map.to.MSA(getLayout(IlluminaBeadChip96Plate), MSA4.plate)


## map a gExperimentSetup to MSA plate see gSetup-class.R




## ## TODO: create generic function
## map.to.MSA <- function(data){
##   names <- colnames(data)
##   nm <- names[names %in% c("chips", "wells")]
##   sampleNames <-  names[!names %in% c("OrigRow","sFactor", "oFactor", "stage", "plates", "chips", "wells", "chipRpws", "chipColumns", "chipRows", "rows", "columns", "chipID", "rowID", "wellID", "cFactor")]
##   if(length(nm) <2 ){stop("need both chips and wells columns in the data")}

##   p <- map.MSA96.BeadChip("IlluminaBeadChip96Plate", MSA4plate="MSA4.plate")
##   p <- p[, c("MSA_ID", nm)]
##   t <- merge(data, p, by=nm, all.x=TRUE)
##   t <- t[order(t$plates, t$chips, t$wells),][,c("plates", "MSA_ID", sampleNames)]
##   return(t)
## }

#t <- map.to.MSA(setup1.data)


# give a list of predefined containers
predefined <- function(){
  l <- c("IlluminaBeadChip", "GenotypingChip",
         "IlluminaBeadChip24Plate",
         "IlluminaBeadChip48Plate",
         "IlluminaBeadChip96Plate",
         "MSA4.plate")
  message("\nAvailable predefined chips and plates.\n\n")
  print(l)
  message("\nType the name of the object for detailed layout information.\n")
}
