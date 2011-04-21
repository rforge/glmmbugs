\name{summaryChain}
\alias{summaryChain}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Compute mean, standard deviation, and quantiles of an MCMC run}
\description{
Computes summary statistics for each parameter.
}
\usage{
summaryChain(chain, probs = c(0.005, 0.025, 0.05, 0.5))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{chain}{The result from  \code{\link{restoreParams}}, or the sims.array component of
   a \code{\link[R2WinBUGS:bugs]{bugs}} call.}
  \item{probs}{Quantiles for the posterior credible interval}
}
\details{
  
}
\value{
 A list of matrices, with rows corresponding to summary statistics and columns to parameters. 
  \item{scalar}{Matrix for the scalar parameters}
  \item{... }{One matrix for each vector valued parameter}
  \item{FittedRateReffect}{For spatial models only, summaries on the natural scale (exponential of FittedReffect).}
}
\references{}
\author{Patrick Brown}
\note{ 
}
\seealso{ \code{\link{restoreParams}} }
\examples{
# create a simple chain
thechain = list(beta = array(1, c(10, 3,4), dimnames = list(NULL, NULL, paste("beta[", 1:4, "]", sep=""))), intercept = matrix(1, 10, 3))

summaryChain(thechain)
}
