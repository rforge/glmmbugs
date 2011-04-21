\name{rongelapResult}
\alias{rongelapResult}
\docType{data}
\title{Rongelap Winbugs Results}
\description{                                  
Results from running Winbugs on the Rongelap data
}
\usage{data(rongelapResult)}
\format{
A list, as produced by the \code{\link[R2WinBUGS]{bugs}} function.
}
\examples{
data(rongelapResult)

rongelapParams = restoreParams(rongelapResult)

}
\seealso{\code{\link{glmmBUGS}}}
\keyword{datasets}
