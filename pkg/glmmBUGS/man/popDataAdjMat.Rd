\name{popDataAdjMat}
\alias{popDataAdjMat}
\docType{data}
\title{Data set containing an adjacency matrix }
\description{
  The \code{popDataAdjMat} Data set contains the adjacency matrix which calculated from the \code{\link[spdep]{poly2nb}} function.
}
\usage{data('popDataAdjMat')}
 
\details{
 It is a adjacency matrix denoting the neighbours of Ontario census subdivisions. 
}
 
\examples{
\dontrun{
library('diseasemapping')
data('popdata')
popDataAdjMat = spdep::poly2nb(popdata,
	row.names=as.character(popdata$CSDUID))
}
data('popDataAdjMat')
summary(popDataAdjMat)
attributes(popDataAdjMat)$region.id[1:10]
}
\keyword{datasets}
