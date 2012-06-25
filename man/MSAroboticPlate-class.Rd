\name{MSAroboticPlate-class}
%\Rdversion{1.1}
\docType{class}
\alias{MSAroboticPlate-class}
\alias{initialize,MSAroboticPlate-method}
%\alias{map.to.MSA,data.frame,MSAroboticPlate-method}
%\alias{map.to.MSA,gExperimentSetup,MSAroboticPlate-method}
\alias{show,MSAroboticPlate-method}
\alias{MSA4.plate}
\alias{BeadChip96ToMSA4MAP}
\alias{getLayout,MSAroboticPlate-method}
\alias{map.to.MSA,data.frame,MSAroboticPlate-method}
\alias{map.to.MSA,gExperimentSetup,MSAroboticPlate-method}
\alias{map.to.MSA}

\title{Class \code{"MSAroboticPlate"}}
\description{
  A class store layout information of MSA robotic loader plate.
  % Run \code{RShowDoc("Meth450_Tracking_Sheet_onepage.pdf",
  % package="OSAT")} for more details.
}
% \section{Objects from the Class}{
% Objects can be created by calls of the form \code{new("MSAroboticPlate", nRows, nColumns, byrow, comment, ...)}.
% %%  ~~ describe objects here ~~ 
% }
% \section{Slots}{
%   \describe{
%     \item{\code{wells}:}{Object of class \code{"character"} ~~ }
%     \item{\code{nRows}:}{Object of class \code{"integer"} ~~ }
%     \item{\code{nColumns}:}{Object of class \code{"integer"} ~~ }
%     \item{\code{layout}:}{Object of class \code{"data.frame"} ~~ }
%     \item{\code{byrow}:}{Object of class \code{"logical"} ~~ }
%     \item{\code{metadata}:}{Object of class \code{"list"} ~~ }
%   }
% }
\section{Extends}{
Class \code{"\linkS4class{gArray}"}, directly.
}
\section{Methods}{
  \describe{
%    \item{initialize}{\code{signature(.Object = "MSAroboticPlate")}: ... }
    \item{map.to.MSA}{\code{signature(x = "data.frame", y =
	"MSAroboticPlate")}: Return a data frame that assign samples
      (represented by the input data frame) to MSA robotic plate sequentially.  }
    \item{map.to.MSA}{\code{signature(x = "gExperimentSetup", y =
	"MSAroboticPlate")}: Return a data frame  that assign samples from an
      \code{gExperimentSetup} object to MSA robotic plate. }
    \item{show}{\code{signature(object = "MSAroboticPlate")}: Shows the
      layout of the plate. }
	 }
}

\section{Predefined objects}{
  
\code{MSA4.plate}: A predeined object of class \code{MSAroboticPlate}
that represent a 96 position plate. 

%data(MSA4.plate)

\code{BeadChip96ToMSA4MAPMap}: The loading order a MSA4 robotic loader
used to load BeadChips.
  }
% \references{
% %%  ~~put references to the literature/web site here~~
% }
% \author{
% %%  ~~who you are~~
% }
% \note{
% %%  ~~further notes~~
% }

% %% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%%  ~~objects to See Also as \code{\link{~~fun~~}}, ~~~
\code{\linkS4class{gExperimentSetup}}
}
\examples{
library("OSAT")
# data as an example
inPath <- system.file("extdata", package="OSAT")
pheno <- read.table(file.path(inPath, 'samples.txt'), header=TRUE, sep="\t")

## create object to hold sample information
gs <- setup.sample(pheno,  optimal=c("SampleType", "Race", "AgeGrp"), strata=c("SampleType") )
gs

gc <- setup.container(IlluminaBeadChip96Plate, 6, batch='plates')
gc

gSetup <- create.optimized.setup(sample=gs, container=gc, nSim=100)
out <- map.to.MSA(gSetup, MSA4.plate)
}
\keyword{classes}
