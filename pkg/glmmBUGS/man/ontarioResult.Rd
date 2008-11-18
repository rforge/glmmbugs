\name{ontarioResult}
\alias{ontarioResult}
\docType{data}
\title{Ontario Winbugs Results}
\description{                                  
Results from running Winbugs on the ontario data
}
\usage{data(ontarioResult)}
\format{
A list, as produced by the \code{\link[R2WinBUGS]{bugs}} function.
}
\examples{
data(ontarioResult)

ontarioParams = restoreParams(ontarioResult)

ontarioSummary = summaryChain(ontarioParams)
}
\keyword{datasets}
