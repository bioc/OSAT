\name{QC}
\alias{QC}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  QC
%%  ~~function to do ... ~~
}
\description{
Provide a visual summary of the sample placement per plates, and
Chi-squre test of dependence between plates and other considered variables
from sample. 
}
\usage{
QC(object, main = NULL, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
    An object of class \code{gExperimentSetup}.
}
\item{main}{
  Mail title on the bar plot.
}
\item{\dots}{
  Additional plot parameters.
}
}
% \details{
% %%  ~~ If necessary, more details than the description above ~~
% }
% \value{
% %%  ~Describe the value returned
% %%  If it is a LIST, use
% %%  \item{comp1 }{Description of 'comp1'}
% %%  \item{comp2 }{Description of 'comp2'}
% %% ...
% }
% \references{
% %% ~put references to the literature/web site here ~
% }
% \author{
% %%  ~~who you are~~
% }
% \note{
% %%  ~~further notes~~
% }

% %% ~Make other sections like Warning with \section{Warning }{....} ~

% \seealso{
% %% ~~objects to See Also as \code{\link{help}}, ~~~
% }
\examples{
library("OSAT")
inPath <- system.file("extdata", package="OSAT")
pheno <- read.table(file.path(inPath, 'samples.txt'), header=TRUE, sep="\t")

## create object to hold sample information
gs <- setup.sample(pheno,  optimal=c("SampleType", "Race", "AgeGrp"), strata=c("SampleType") )
## create object that represents the used in the experiment
gc <- setup.container(IlluminaBeadChip96Plate, 6, batch='plates')
## create an optimized setup. 
# demonstration only. nSim=5000 or more are commonly used.
gSetup <- create.optimized.setup(sample=gs, container=gc, nSim=500)
QC(gSetup)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
