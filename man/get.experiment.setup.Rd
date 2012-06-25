\name{get.experiment.setup}
\alias{get.experiment.setup}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Get experiment assignment after optimization.
}
\description{
Return a data frame that contains sample and assigned well postion.
}
\usage{
get.experiment.setup(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ A \code{gExperimentSetup} object. }
}
\value{
A data frame is returned in the order of initial sample data frame, with
additional columns indicate assgined well position.
}
% \references{
% %% ~put references to the literature/web site here ~
% }
% \author{
% %%  ~~who you are~~
% }
% \note{
% %%  ~~further notes~~
% }

%% ~Make other sections like Warning with \section{Warning }{....} ~

% \seealso{

% }
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
# demonstration only. nSim=5000 or more are commonly used.
gSetup <- create.optimized.setup(sample=gs, container=gc, nSim=500)

mySetup <- get.experiment.setup(gSetup)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
