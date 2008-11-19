\name{ontario}
\alias{ontario}
\docType{data}
\title{Ontario data on molar cancer}
\description{
Data frame showing expected and observed counts of molar cancer in Ontario
}
\usage{data(ontario)}
\format{
  A data frame with 585 observations on the following 3 variables.
  \describe{
    \item{\code{CSDUID}}{factor of Ontario census subdivision ID numbers}
    \item{\code{Observed}}{Observed molar cancer cases}
    \item{\code{logExpected}}{expected cases}
  }
}
\details{
See the documentation for \link{glmmBUGS} for how this was created. 
}
\examples{
data(ontario)
head(ontario)
}
\keyword{datasets}
