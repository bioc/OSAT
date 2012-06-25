\name{optimal.shuffle}
\alias{optimal.shuffle}

\title{
  optimal.shuffle
%%   ~~ data name/kind ... ~~
}
\description{
Optimize a sample assingment setup by shuffling.
}
\details{
  Given any \code{gExperimentSetup} object, we randomly select k samples from different batches and shuffle them between batches to create a new sample assignment. k = 2 by default but could be any number up to half of the sample size. Value of the objective function is calculated on the new setup and compared to that of the original one. If the value is smaller then the new assignment replaces the previous one. This procedure will continue until we reach a preset number of attempts (usually in the tens of thousands).
  }
\usage{optimal.shuffle(x, nSim, k)}
\arguments{
  \item{x}{A \code{gExperimentSetup} object.
}
  \item{nSim}{Number of shuffling steps.
  }
  \item{k}{Number of samples been shuffled. Default \code{k=2} when omitted..
  }
}
\value{
A class \code{gExperimentSetup} object, after optimized.
}
\seealso{
  \code{\link{optimal.block}}
%% ~~objects to See Also as \code{\link{help}}, ~~~
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
# demonstration only. nSim=5000 or more are commonly used.
g3 <- optimal.shuffle(gSetup0, nSim=500, k=2)
}
% \keyword{function}
