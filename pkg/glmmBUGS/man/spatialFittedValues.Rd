\name{spatialFittedValues}
\alias{spatialFittedValues}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Compute summaries of fitted values for geostatistical models
}
\description{
Add intercept to spatial random effects and exponentiate, compute summaries of posterior.
}
\usage{
spatialFittedValues(chain, threashold=1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{chain}{
output of restoreParams function
}
\item{threashold}{value which the posterior probability of exceeding is computed.  Defaults to 1}
}
\details{
}
\value{A list with posterior mean, posterior standard deviation, and posterior probability of
exceeding 1 for exp(mu + U(s))
}
\references{
}
\author{
Patrick Brown
}
\note{
}


\seealso{
}
\examples{
\dontrun{
    rongelapSummary$fitted = spatialFittedValues(rongelapParams)

  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
