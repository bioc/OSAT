\name{create.optimized.setup}
\alias{create.optimized.setup}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Create optimized setup
}
\description{
Create a optimized sample assignment.
}
\usage{
create.optimized.setup(fun = "default", sample, container, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fun}{ The name of the optimization function. When omitted, it is the same
   as \code{default} or \code{optimal.shuffle}. Currently
   alternative optimization function is \code{optimal.block}. User can
   define their own optimization function.
}
  \item{sample}{ A \code{gSample} object contains sample related information.
   
}
  \item{container}{ A \code{gContainer} object contains experiment
  related container (chips, plate, exclusion etc) information.
}
\item{\dots}{ Additional parameters passed to optimization function.
}
}
\details{
 Currently two methods are
available: \code{optimal.shuffle} (default) and \code{optimal.block}
methods. The function is equivalent to \code{create.experiment.setup()}
followed by corresponding optimization function of the same name.
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
A \code{gExperimentSetup} object is returned to store all related information.
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
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
