\name{gContainer-class}

\docType{class}
\alias{gContainer-class}
\alias{cid,gContainer-method}
\alias{get.gAssembly,gContainer-method}
\alias{getLayout,gContainer-method}
\alias{initialize,gContainer-method}
\alias{metadata<-,gContainer-method}
\alias{metadata,gContainer-method}
\alias{show,gContainer-method}
\alias{setup.container}
\alias{exclude<-}
\alias{exclude<-,gContainer-method}

\title{Class \code{"gContainer"}}
\description{
A class for storage information related to the experiment container
setup. 
}
\section{Objects from the Class}{
  Objects can be created by calls of the form
  \code{setup.container(plate, n, batch, exclude)}.

}
\usage{
setup.container(plate, n, batch = "plates", exclude = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{plate}{A object of \code{gPlate} class or its expanded
    class. 
}
  \item{n}{Number of plates used in the experiment.
}
  \item{batch}{The level where batch effect is considered.
}
  \item{exclude}{ A data frame indicate location of wells in the
    container that should be  excluded from sample
    assigment in the experiment.
}

}
% \details{
% %%  ~~ If necessary, more details than the description above ~~
% }
\value{
A \code{gContainer} object.
}       
\section{Slots}{
  \describe{
    \item{\code{plate}:}{Object of class \code{"gPlate"}. Plate used in
      the experiment.  }
    \item{\code{n}:}{Object of class \code{"integer"}. Number of plate used. }
    \item{\code{batch}:}{Object of class \code{"character"}. On what
      level batch effect are considered. Could be \code{"plates"} or
      \code{"chips"}. }
    \item{\code{exclude}:}{Object of class \code{"data.frame"}. A data
      frame indicate wells that should be excluded from the sample Assignment. }
    \item{\code{data}:}{Object of class \code{"list"}. A list holds
      summaries and other useful informaiton.  }
    \item{\code{metadata}:}{Object of class \code{"list"}. A list for a brief
      description and other useful informtion.}
  }
}
\section{Methods}{
  \describe{
%     \item{cid}{\code{signature(x = "gContainer")}: ... }
    \item{get.gAssembly}{\code{signature(x = "gContainer")}: Return the
      assembly of the plates used in the experiments. }
    \item{getLayout}{\code{signature(x = "gContainer")}:  Return a data
      frame that holds the layout of available wells in the container. }
    %     \item{initialize}{\code{signature(.Object = "gContainer")}: ... }
    \item{exclude<-}{\code{signature(object = "gContainer")}: Exclude some
      wells/chips/plates from the container for randomization. }
    \item{metadata<-}{\code{signature(x = "gContainer")}: Set the metadata. }
    \item{metadata}{\code{signature(x = "gContainer")}: Get the metadata. }
    \item{show}{\code{signature(object = "gContainer")}: A brief summary
    of the container.}
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
library("OSAT")
# a container consist of 6 predefined Illumina plates with 96 wells each
gc <- setup.container(IlluminaBeadChip96Plate, 6, batch='plates')
gc

# to exclude first wells on first chips of each plate.
(excludedWells <- data.frame(plates=1:6, chips=rep(1,6), wells=rep(1,6)) )
gc3 <- setup.container(IlluminaBeadChip96Plate, 6, batch='plates', exclude=excludedWells)
}
\keyword{classes}
