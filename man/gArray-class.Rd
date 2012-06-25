\name{gArray-class}

\docType{class}
\alias{gArray-class}
\alias{dim,gArray-method}
\alias{getLayout,gArray-method}
\alias{metadata<-,gArray-method}
\alias{metadata,gArray-method}
\alias{ncol,gArray-method}
\alias{nrow,gArray-method}
\alias{getLayout}
\alias{metadata<-}
\alias{metadata}

\title{Class \code{"gArray"}}
\description{
  The virtual class \code{gArray} is used to stadarize the storage and  access of array-like experiment arrangement.It represent the  arrangement of the experiment units, such as slide ans plate.
%%  ~~ A concise (1-5 lines) description of what the class is. ~~
}
\section{Objects from the Class}{A virtual Class: No objects may be created from it.}
\section{Slots}{
  \describe{
    \item{\code{nRows}:}{Object of class \code{"integer"}. Number of
      units per row. }
    \item{\code{nColumns}:}{Object of class \code{"integer"}. Number of
      units per column. }
    \item{\code{layout}:}{Object of class \code{"data.frame"}. A data
      frame represent the layout and order of postiion. }
    \item{\code{byrow}:}{Object of class \code{"logical"}. The order of
      position is filled by row or not. Default is FALSE, same as
      used in mactrix().  }
    \item{\code{metadata}:}{Object of class \code{"list"}. A list that
      holds brief description of the class and other information. }
  }
}
\section{Methods}{
  In the following code snippets, \code{x} is an \code{gArray} object.
  
  \describe{
    \item{dim}{\code{signature(x = "gArray")}: Return the dimmension of
    the array. }
    \item{getLayout}{\code{signature(x = "gArray")}: Return a data frame
    represent the layout postions of the array.}
    \item{metadata<-}{\code{signature(x = "gArray")}: Set the metadata. }
    \item{metadata}{\code{signature(x = "gArray")}: Get the metadata. }
    \item{ncol}{\code{signature(x = "gArray")}: Return the number of columns.}
    \item{nrow}{\code{signature(x = "gArray")}: Return the number of rows. }
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
% \examples{
% showClass("gArray")
% }
\keyword{classes}
