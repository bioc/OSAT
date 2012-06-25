\name{gSlide-class}
\docType{class}
\alias{gSlide-class}
\alias{getLayout,gSlide-method}
\alias{initialize,gSlide-method}
\alias{show,gSlide-method}

\title{Class \code{"gSlide"}}
\description{
A class repressent a slide used in experiment. Inherited from
Class \code{"gArray"}. The layout of the slide is descriped by a data
frame with columns of "rows", "columns", and "wells".
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("gSlide", nRows, nColumns, byrow, comment)}.
%%  ~~ describe objects here ~~ 
}
%% \section{Slots}{
%%    \describe{
%%      \item{\code{nRows}:}{Object of class \code{"integer"} ~~ }
%%      \item{\code{nColumns}:}{Object of class \code{"integer"} ~~ }
%%      \item{\code{layout}:}{Object of class \code{"data.frame"} ~~ }
%%      \item{\code{byrow}:}{Object of class \code{"logical"} ~~ }
%%      \item{\code{comment}:}{Object of class \code{"character"} ~~ }
%%    }
%% }
\section{Extends}{
Class \code{"\linkS4class{gArray}"}, directly.
}
\section{Methods}{
  \describe{
    \item{getLayout}{\code{signature(x = "gSlide")}: ... }
    \item{initialize}{\code{signature(.Object = "gSlide")}: ... }
    \item{show}{\code{signature(object = "gSlide")}: ... }
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
newSlide <- new("gSlide",  nRows=2, nColumns=2, byrow=FALSE,comment="mock slide")
newSlide
}
\keyword{classes}
