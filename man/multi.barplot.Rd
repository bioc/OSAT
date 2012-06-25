\name{multi.barplot}
\alias{multi.barplot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  multi.barplot
%%  ~~function to do ... ~~
}
\description{
Plot multiple bar plots based on a single data frame.
}
\usage{
multi.barplot(x, grpVar = "plates", varList, main = NULL, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ A data frame.
}
  \item{grpVar}{ The variable for x-axix.
}
  \item{varList}{A vector of variables, each will be used as y-axis for
   the multiple bar plot.
}
  \item{main}{The overall title of the plot.
}
  \item{\dots}{Additinal parameters for the plot function.
}
}
% \details{
% %%  ~~ If necessary, more details than the description above ~~
% }
\value{
Use the side effect for plot. No return value.
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

% %% ~Make other sections like Warning with \section{Warning }{....} ~

% \seealso{
% %% ~~objects to See Also as \code{\link{help}}, ~~~
% }
\examples{
## create a random assignment and check it
library("OSAT")
data(example.setup)
set.seed(10)
c1 <- getLayout(gc)                 # available wells
c1 <- c1[order(runif(nrow(c1))),]   # shuffle randomly
ranomSetup <- cbind(pheno, c1[1:nrow(pheno), ]) # create a sample assignment

multi.barplot(ranomSetup, grpVar='plates', varList=c("SampleType", "Race", "AgeGrp"), main="A random case")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
