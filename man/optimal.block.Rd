\name{optimal.block}
\alias{optimal.block}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  optimal.block
%%  ~~function to do ... ~~
}
\description{
Optimize a sample assingment setup by selecting from multiple candidate setup.
}
\usage{
optimal.block(x, nSim = 100)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{A \code{gExperimentSetup} object.
}
  \item{nSim}{Number of candidate setup created.
}
}
\details{
Multiple (typically thousands of or more) sample assignment setups are first generated, based only on the list of specified blocking variable(s). Then, the optimal setup is chosen by selecting the setup of sample assignment (from the pool generated in blocking step) which minimizes the value of the objective function based on all variables considered. 
}
\value{
A \code{gExperimentSetup} object, after optimization.
}

\seealso{
  \code{\link{optimal.shuffle}}
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

g2 <- optimal.block(gSetup0, nSim=100)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
