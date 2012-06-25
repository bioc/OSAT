\name{multi.chisq.test}
\alias{multi.chisq.test}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  multi.chisq.test
%%  ~~function to do ... ~~
}
\description{ Do a few Chi-squre tests based on  the same  data frame.
}
\usage{
multi.chisq.test(x, grpVar = "plates", varList, main = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{A data frame.
}
  \item{grpVar}{Common variables. Default is 'plate'.
}
  \item{varList}{A vector of variables. 
}
  \item{main}{The overall title.
}
}
% \details{
% %%  ~~ If necessary, more details than the description above ~~
% }
% \value{
% %%  ~Describe the value returned
% %%  If it is a LIST, use
% %%  \item{comp1 }{Description of 'comp1'}
% %%  \item{comp2 }{Description of 'comp2'}
% %% ...
% }
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

multi.chisq.test(ranomSetup, grpVar='plates', varList=c("SampleType", "Race", "AgeGrp"))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
