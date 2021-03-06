\name{checkChain}
\alias{checkChain}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Plot an MCMC run}
\description{
 Makes time series plots of the parameters (not the random effects) of an MCMC run.
}
\usage{
checkChain(chain, parameters=NULL, oneFigure=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{chain}{The result from  \code{\link{restoreParams}}, or the sims.array componento f a \code{\link[R2WinBUGS:bugs]{bugs}} call.}
\item{parameters}{Vector of character strings giving names of parameters to plot.  Default is all parameters  with names starting with either "beta", "intercept", or "SD".
}
\item{oneFigure}{if \code{TRUE}, use \code{par(mfrow=c(a,b))} to put all plots 
on the same device. Otherwise create a new device for each plot.}
}
\value{
   Plots are produced, nothing is returned
}
 
\author{Patrick Brown}
 
\seealso{ \code{\link{restoreParams}},  \code{\link{summaryChain}}  }
\examples{
% create a sample MCMC result
thechain = list(beta = array(1, c(10, 3,4), 
dimnames = list(NULL, NULL, paste("beta[", 1:4, "]", sep=""))),
 intercept = matrix(1, 10, 3))

checkChain(thechain)

}
