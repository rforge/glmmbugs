\name{rongelapUTM}
\alias{rongelapUTM}
\docType{data}
\title{
Rongelap data 
}
\description{
A SpatialPointsDataFrame containing the Rongelap data, in a UTM projection.
}
\usage{data('rongelapUTM')}
 
\details{
These coordinates were obtained by translating and rotating the original Rongelap data
until all the coordinates fit into the Rongelap border given by www.gadm.org.  
So they are not exact.
}
\source{
See the help file for \code{\link[geoRglm]{rongelap}, or 
\url{http://www.leg.ufpr.br/doku.php/pessoais:paulojus:mbgbook:datasets}
}
\examples{
data('rongelapUTM')
plot(rongelapUTM)
\dontrun{
library('raster')
rongelapBorderLL = getData('GADM', country='MHL',level=0) 
library('rgdal')
rongelapBorderUTM = spTransform(rongelapBorderLL, 
	CRS(proj4string(rongelapUTM)))
plot(rongelapBorderUTM, add=TRUE)
}
}
\keyword{datasets}
