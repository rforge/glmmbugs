\name{rongelapUTM}
\alias{rongelapUTM}
\docType{data}
\title{
Rongelap data 
}
\description{
A SpatialPointsDataFrame containing the Rongelap data, in a UTM projection.
}
\usage{data(rongelapUTM)}
\format{
}
\details{
These coordinates were obtained by translating and rotating the original Rongelap data
until all the coordinates fit into the Rongelap border given by www.gadm.org.  
So they are not exact.
}
\source{
See the help file for rongelap in geoRglm, or http://www.leg.ufpr.br/doku.php/pessoais:paulojus:mbgbook:datasets
}
\references{
}
\examples{
data(rongelapUTM)
\dontrun{
library(rgdal)
rongelapLL<-spTransform(rongelapUTM,CRS("+proj=longlat +datum=NAD83"))
   load(url("http://www.gadm.org/data/rda/MHL_adm0.RData"))
   plot(rongelapLL)
   plot(gadm, add=T)    
}
}
\keyword{datasets}
