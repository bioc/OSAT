\name{gAssembly-class}
\docType{class}
\alias{gAssembly-class}
\alias{getLayout,gAssembly-method}
\alias{initialize,gAssembly-method}
\alias{metadata,gAssembly-method}
\alias{metadata<-,gAssembly-method}
\alias{show,gAssembly-method}

\title{Class \code{"gAssembly"}}
\description{
A class represent an assembly of experiment plates.
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("gAssembly", plate, n, comment)}.
%%  ~~ describe objects here ~~ 
}
\section{Slots}{
  \describe{
    \item{\code{plate}:}{Object of class \code{"gArray"}. Plate used in
      the experiment.  }
    \item{\code{n}:}{Object of class \code{"integer"}. NUmber of plates used. }
    \item{\code{layout}:}{Object of class \code{"data.frame"}. The
      layout of wells in the assembly.}
    \item{\code{metadata}:}{Object of class \code{"list"}. A list holds
      meta data of the assembly.}
  }
}

\section{Methods}{
  \describe{
    \item{getLayout}{\code{signature(x = "gAssembly")}: Return a data
  frame that represent the wells location in the assembly. }
%    \item{initialize}{\code{signature(.Object = "gAssembly")}: }
\item{metadata}{\code{signature(x = "gAssembly")}: Get metadata. }
\item{metadata<-}{\code{signature(x = "gAssembly")}: Set metadata. }
    \item{show}{\code{signature(object = "gAssembly")}: Print a brief
  summary of the assembly. }
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

\seealso{
  \code{"\linkS4class{gContainer}"}.
}
% \examples{
% showClass("gAssembly")
% }
\keyword{classes}
