## genetic functions used in the OSAT package
## Auther: Li Yan
## 

if (!isGeneric("metadata"))
  setGeneric("metadata", function(x, ...) standardGeneric("metadata"))

if (!isGeneric("metadata<-"))
setGeneric("metadata<-",
           function(x, ..., value) standardGeneric("metadata<-"))

if (!isGeneric("nrow"))
    setGeneric("nrow", function(x) standardGeneric("nrow"))
if (!isGeneric("ncol"))
    setGeneric("ncol", function(x) standardGeneric("ncol"))
## dim is S3 method
##if (!isGeneric("dim")) setGeneric("dim")
#    setGeneric("dim", function(x) standardGeneric("dim"))

## layout is a common function in graphic packages, so use diff names
if (!isGeneric("getLayout"))
    setGeneric("getLayout", function(x, ...) standardGeneric("getLayout"))

if (!isGeneric("cid"))
    setGeneric("cid", function(x) standardGeneric("cid"))

if (!isGeneric("get.gAssembly"))
    setGeneric("get.gAssembly", function(x) standardGeneric("get.gAssembly"))

## map a list of positions in the given order to the MSA plate
setGeneric("map.to.MSA", function(x, y, ...) standardGeneric("map.to.MSA"))

## retrive samples 
#if (!isGeneric("samples"))
    setGeneric("samples", function(x) standardGeneric("samples"))

if (!isGeneric("sid"))
    setGeneric("sid", function(x) standardGeneric("sid"))

setGeneric("summary")  ## this will turn S3 generic to S4 generic


## # define a generic function for object function?
## setGeneric("get.experiment.setup", function(x, ...) standardGeneric("get.experiment.setup"))
## setMethod("get.experiment.setup", "gContainer", function(x, ...) x@layout)



if (!isGeneric("plot"))
      setGeneric("plot", function(x, y, ...) standardGeneric("plot"))


## setGeneric("setup.sample", function(x, optimal, strata, ...) standardGeneric("setup.sample"))

setGeneric("exclude<-",
           function(x,value) standardGeneric("exclude<-"))

