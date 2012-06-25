\name{gExperimentSetup-class}

\docType{class}
\alias{gExperimentSetup-class}
\alias{get.gAssembly,gExperimentSetup-method}
\alias{initialize,gExperimentSetup-method}

%\alias{plot,gExperimentSetup,missing-method}
\alias{samples,gExperimentSetup-method}
\alias{show,gExperimentSetup-method}
\alias{summary,gExperimentSetup-method}
\alias{create.experiment.setup}
\alias{metadata,gExperimentSetup-method}
\alias{metadata<-,gExperimentSetup-method}
%\alias{map.to.MSA}

\title{Class \code{"gExperimentSetup"}}
\description{
A class stores all relevant informations related to the experimantal sample assgnment.
}
\details{
Function \code{create.experiment.setup()} create a block randomized
experiment assignment, without optimization step.

Functions \code{optimal.shuffle()} and \code{optimal.block} optimize
setup using different optimization methods.

Functiion \code{create.optimized.setup()} create a optimized setup directly.
}

\section{Slots}{
  \describe{
    \item{\code{expSetup}:}{Object of class \code{"data.frame"}. A data
  frame represents the sample placement to well locations in the container. }
    \item{\code{data}:}{Object of class \code{"list"}.  A list that hold
  a \code{gSample} object, a \code{gContainer} object and assignment
  link the two objects.}
    \item{\code{summaryInfo}:}{Object of class \code{"list"}. Some
  summary of the sample, container and assignment. }
    \item{\code{metadata}:}{Object of class \code{"list"}. Metadata of
  the object. }
  }
}
\section{Accessors}{
  \describe{
    \item{get.gAssembly}{\code{signature(x = "gExperimentSetup")}: Get the chip/plate assembly used for container.}
%    \item{initialize}{\code{signature(.Object = "gExperimentSetup")}: ... }
    \item{samples}{\code{signature(x = "gExperimentSetup")}: Return the
      \code{gSample} object. }
    \item{get.gContainer}{: ...}
    \item{metadata}{\code{signature(x = "gExperimentSetup")}: Get the metadata. }
    \item{metadata<-}{\code{signature(x = "gExperimentSetup")}: Set the metadata. }
    
  }
}
\section{Related Methods}{
  \describe{
    \item{map.to.MSA}{\code{signature(x = "gExperimentSetup", y =
  "MSAroboticPlate")}: A method map the experiment setup to a set of MSA 96
  wells robotic plates. }
    \item{plot}{\code{signature(x = "gExperimentSetup", y = "missing")}:
  Visual  }
    \item{show}{\code{signature(object = "gExperimentSetup")}: ... }
    \item{summary}{\code{signature(object = "gExperimentSetup")}: ... }
    \item{get.experiment.setup}{: Return a data frame with linked sample
  and container information.}
	 }
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

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%%  ~~objects to See Also as \code{\link{~~fun~~}}, ~~~
\code{\linkS4class{MSAroboticPlate}} 
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

gSetup0 <- create.experiment.setup(sample=gs, container=gc)

}
\keyword{classes}
