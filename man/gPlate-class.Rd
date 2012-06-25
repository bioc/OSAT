\name{gPlate-class}

\docType{class}
\alias{gPlate-class}
\alias{initialize,gPlate-method}
\alias{show,gPlate-method}
\alias{getLayout, gPlate-method}

\title{Class \code{"gPlate"}}
\description{
A virtual class represent a experiment plate that is consistent of multiple slides.
}
\section{Objects from the Class}{A virtual Class: No objects may be created from it.}
\section{Slots}{
  \describe{
    \item{\code{chip}:}{Object of class \code{"gSlide"}. A virtual class represent a experiment plate that is consistent of multiple slides. }
    \item{\code{nRows}:}{Object of class \code{"integer"}. Number of
      slide per row.  }
    \item{\code{nColumns}:}{Object of class \code{"integer"}. Number of
      slide per column.  }
    \item{\code{layout}:}{Object of class \code{"data.frame"}. A data
      frame represent the layout of the plate. Usually down to the wells
    level. }
    \item{\code{byrow}:}{Object of class \code{"logical"}. The order of
      the position is filled by row or not. Default is FALSE, same as
      used in mactrix().  }
    \item{\code{metadata}:}{Object of class \code{"list"}. A list that
      holds brief description of the class and other information.   }
  }
}
\section{Extends}{
Class \code{"\linkS4class{gArray}"}, directly.
}
% \section{Methods}{
%   \describe{
%     \item{initialize}{\code{signature(.Object = "gPlate")}: ... }
%     \item{show}{\code{signature(object = "gPlate")}: ... }
% 	 }
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
% \examples{
% showClass("gPlate")
% }
\keyword{classes}
