\name{gSample-class}
\docType{class}
\alias{gSample-class}
\alias{show,gSample-method}
\alias{setup.sample}
\alias{summary}
\alias{summary,gSample-method}
\alias{initialize,gSample-method}

\title{Class \code{gSample}}
\description{
Create a class used for storage of sample related information.
}
% \section{Objects from the Class}{
% Objects can be created by calls of the form \code{new("gSample", ...)}.
% %%  ~~ describe objects here ~~ 
% }
\usage{
setup.sample(x, optimal, strata)
summary(object, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ a data frame holds sample variables.
}
  \item{optimal}{ a vector of sample varible names to be treated as optimal variables.
}
  \item{strata}{ a vector of sample varible names to be treated as block
    variables. If omitted, the first element of \code{optimal} vector is
    treated as strata.
  }
  \item{object}{ an object of class \code{gSample}.}
  \item{...}{additional arguments affecting the summary produced.}
}
% \details{
% %%  ~~ If necessary, more details than the description above ~~
% }
\value{
  An object of \code{gSample} class for function \code{setup.sample()}.
  A list of two tables for function \code{summary(object)}:
  \item{strataTable}{Sample freqency table by stata variables.}
  \item{optimalTable}{Sample freqency table by optimal variables.}
} 



\section{Slots}{
  \describe{
    \item{\code{rawData}:}{Object of class \code{"data.frame"} The
      original data frame that holds the sample informaiton. }
    \item{\code{optimal}:}{Object of class \code{"character"} A
      charactor vector of column names in the sample info data
      frame. Represent the optimal variables. }
    \item{\code{strata}:}{Object of class \code{"character"}  A
      charactor vector of column names in the sample info data
      frame. Represent the blocking variables. }
    \item{\code{data}:}{Object of class \code{"list"} A list holds
      summaries and other useful information.}
  }
}
\section{Methods}{
  \describe{
    \item{show}{\code{signature(object = "gSample")}: ... }
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

% %% ~Make other sections like Warning with \section{Warning }{....} ~

% \seealso{
% %%  ~~objects to See Also as \code{\link{~~fun~~}}, ~~~
% %%  ~~or \code{\linkS4class{CLASSNAME}} for links to other classes ~~~
% }
\examples{
inPath <- system.file("extdata", package="OSAT")
pheno <- read.table(file.path(inPath, 'samples.txt'), header=TRUE, sep="\t")
## create object to hold sample information
gs <- setup.sample(pheno,  optimal=c("SampleType", "Race", "AgeGrp"), strata=c("SampleType") )
}
\keyword{classes}
