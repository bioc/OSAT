\name{BeadChip-class}

\docType{class}
\alias{BeadChip-class}
\alias{IlluminaBeadChip}
\alias{GenotypingChip}
\alias{getLayout,BeadChip-method}

\title{Class \code{"BeadChip"}}
\description{
  A predefined class represent BeadChip from Illumina. Inherited from
  class \code{gSlide}.
  %% \textregistered.
}
%% \section{Objects from the Class}{
%% Objects can be created by calls of the form \code{new("BeadChip", nRows, nColumns, byrow, comment)}.
%%  ~~ describe objects here ~~ 
%% }
\section{Predefined chips}{
  \code{IlluminaBeadChip}: A slide with 6 rows and 2 columns. The 12
  wells is numbered columns first.
  
  \code{GenotypingChip}: A slide with 12 rows and one column.
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
% \usage{
% newChip <- new("BeadChip",  nRows=6, nColumns=1, byrow=FALSE, comment="mock chip")
% data(GenotypingChip)
% data(IlluminaBeadChip)
% }
\section{Extends}{
Class \code{"\linkS4class{gSlide}"}, directly.
Class \code{"\linkS4class{gArray}"}, by class "gSlide", distance 2.
}
\section{Methods}{
  \describe{
    \item{getLayout}{\code{signature(x = "BeadChip")}: Return a data
  frame represent the layout of the chip. }
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
%library("OSAT")
\examples{
newChip <- new("BeadChip",  nRows=6, nColumns=1, byrow=FALSE, comment="mock chip")
newChip

GenotypingChip
IlluminaBeadChip
}

\keyword{classes}
