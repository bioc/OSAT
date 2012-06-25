\name{BeadPlate-class}
\docType{class}
\alias{BeadPlate-class}


\alias{IlluminaBeadChip96Plate}
\alias{IlluminaBeadChip48Plate}
\alias{IlluminaBeadChip24Plate}

\title{Class \code{"BeadPlate"}}
\description{
A class represent a plate consistent of multiple chips. Inherite from
class \code{gPlate}. 
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("BeadPlate", chip, nRows, nColumns, comment, ...)}.
%%  ~~ describe objects here ~~ 
}
\section{Slots}{
  \describe{
    \item{\code{chip}:}{Object of class \code{"gSlide"}.  The slide used
      in the plate.}
    \item{\code{nRows}:}{Object of class \code{"integer"}. Number of
      clide per row. }
    \item{\code{nColumns}:}{Object of class \code{"integer"}. Number of
      slide per column. }
    \item{\code{layout}:}{Object of class \code{"data.frame"}. A data frame represent
          the layout of the plate.   }
    \item{\code{byrow}:}{Object of class \code{"logical"}. Order of
      slide on the plate. }
    \item{\code{metadata}:}{Object of class \code{"list"}.A brief description of the plate.  }
  }
}
\section{Predefined plate}{
  \code{IlluminaBeadChip96Plate}
  
  \code{IlluminaBeadChip48Plate}
  
  \code{IlluminaBeadChip24Plate}
  
  are plates that hold 2, 4, 8 \code{IlluminaBeadChip} chips and have 24, 48, 96 wells, respectively.
  }
\section{Extends}{
Class \code{"\linkS4class{gPlate}"}, directly.
Class \code{"\linkS4class{gArray}"}, by class "gPlate", distance 2.
}
% \section{Methods}{
% No methods defined with class "BeadPlate" in the signature.
% }
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
newPlate <- new("BeadPlate", chip=IlluminaBeadChip, nRows=2L, nColumns=4L,  comment="mock plate")
newPlate
}
\keyword{classes}
