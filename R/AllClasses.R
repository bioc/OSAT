


## define classes representing containers (chip, plate etc)
## used in genomic experiments, and classes representing the assembly
## of all plate/chips/wells
## Auther: Li Yan
## 


## we should NOT derive it from matrix since we do not want to do any mtrix operation on it. However many aspect of matrix can be borrowed.

## VIRTUAL class, anything that have a rectangular layoout
## TODO: validation function?
## byrow the same as in matrix()
setClass("gArray",
         representation(nRows="integer", nColumns="integer",
                        layout="data.frame", byrow="logical",
                        metadata="list",
                        "VIRTUAL"),
          prototype=prototype(nRows=1L, nColumns=1L, byrow=FALSE,
            metadata=list(comment=character())
            )
         )
# for debug resetClass() findClass()


## any Slide can derived from here
## 
setClass("gSlide", contains="gArray",
          prototype=prototype(nRows=6L, nColumns=2L, byrow=FALSE
            )
         )

## Illumina BeadChip class
setClass("BeadChip",
         contains="gSlide",
         prototype=prototype(metadata=list(comment="By Illumina"))
         )



## virtual class for plate 
## setClass("gPlate",
##          representation(chip="gSlide", "VIRTUAL"),
##          contains="gArray",
##            prototype=prototype(chip=IlluminaBeadChip, nRows=2L, nColumns=4L
##              )
##          )
# modify initialization, do not use prototype
setClass("gPlate",
         representation(chip="gSlide", "VIRTUAL"),
         contains="gArray")

## BeadChip Plate class
setClass("BeadPlate",
         contains="gPlate")


## assembly class
## extend an gArray 
## could be plate or chip
## TODO: chip level not working yet
## TODO: more flexible assembly
## define more general container?
## TODO: more metadata? move plate/n/ etc to metadata and move container out?
## setClass("gAssembly",
##          representation(plate="gArray",
##                         n="integer",
##                         layout="data.frame", 
##                         metadata="list"),
##           prototype=prototype(plate=IlluminaBeadChip96Plate, n=2L
##             )
##          )
setClass("gAssembly",
         representation(plate="gArray",
                        n="integer",
                        layout="data.frame", 
                        metadata="list")
         )
## not working yet
# tmp1 <- new("gAssembly", plate=GenotypingChip, n=2, comment="test chip")


## TODO: should gContainer inherite gAssembly?
# data slot include assembly, cid, container, blockTable
#   assembly: assembly created based on plates and number
#   container: layout of assembly removed excluded position, add blocking factor
#   cid: only wellID and cfactor
#   blockTable: freq count of available container position by block
setClass("gContainer",
         representation(plate="gPlate", n="integer",
                        batch="character",
                        exclude="data.frame",
                        data="list",
                        metadata="list"
                        ),
         prototype=prototype(exclude=data.frame(), metadata=list(comment=character()) )
         )


## TODO: use setIs()   ??
## TODO: use callNextMethod


## MSA robatic plate can derived from here
## wells: represent the order and name of the wells
## TODO: BC possition? not sure
setClass("MSAroboticPlate", contains="gArray",
         representation("wells"="character"),
          prototype=prototype(nRows=8L, nColumns=12L, byrow=FALSE
            )
         )


# define a class representing the list of samples in a experiment
# Include:
#         response variable (continous (not yet), categorical)
#         list of categorical independent variables (factors)
#         list of continouse independent variables (optional)
# Method:
#         init / creation
#         basic summaries
#         validation (later)
# only deal with relatively small sample size, < 10000?

setClass("gSample",
         representation=representation(rawData='data.frame',
           optimal='character',
           strata='character',          
           data='list'
           ),
         prototype=prototype(
           data=list(sid=NULL,
             annotation=NULL, strataTable=NULL, optimalTable=NULL)
           )
         )


## define classes representing experiment sample assignment

setClass("gExperimentSetup",
         representation(expSetup="data.frame",
                        data="list",
                        summaryInfo='list',
                        metadata='list'
                        ),
         prototype=prototype(data=list(sample=NULL, container=NULL, link=NULL),
                             summaryInfo=list(average.block=matrix(),
                                              average.optimal=matrix()),
                             metadata=list(optimalFunction=NULL, optValue=NULL)
            )
         )
# resetClass("gExperimentSetup", "OSAT")

