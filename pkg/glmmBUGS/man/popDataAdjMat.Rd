\name{popDataAdjMat}
\alias{popDataAdjMat}
\docType{data}
\title{Data set contains the adjacency matrix }
\description{
  The popDataAdjMat Data set contains the adjacency matrix which calculated from the poly2nb function.
}
\usage{data(popDataAdjMat)}
\format{
}
\details{
 It is a adjacency matrix denoting the neighbours of Ontario census subdivisions. Created by:
 \code{library(diseasemapping); data(popdata); popDataAdjMat = poly2nb(ontario,row.names=as.character(ontario[["CDSUID"]]))
}. 
}
\source{

}
\references{
 
}
\examples{
data(popDataAdjMat)
}
\keyword{datasets}
