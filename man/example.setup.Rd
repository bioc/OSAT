\name{example.setup}
\alias{gSetup}
\alias{gSetup2}
\alias{gs}
\alias{gc}
\alias{pheno}

\docType{data}
\title{
R data object based on example file.
}
\description{Data objects created using the example data for demonstration.}
\usage{data(example.setup)}
\section{Details}{
\describe{
  See vignette for codes that created these data objects.
  \item{pheno:}{ a data frame created from example sample information text file.}
  \item{gs:}{ a class \code{gSample} object created based on example
    sample using function \code{setup.sample().}.}
  \item{gc:}{ a class \code{gContainer} object using function \code{setup.container()}.}
  \item{gSetup:}{ a class \code{gExperimentSetup} object. Created using
    \code{optimal.shuffle()} method.}
  \item{gSetup2:}{ a class \code{gExperimentSetup} object. Created using
    \code{optimal.block()} method.}
}
}
\examples{
library("OSAT")
data(example.setup)
head(pheno)
gs
gc
gSetup
}
\keyword{datasets}
