\name{rongelapUTM}
\alias{rongelapUTM}
\alias{rongelapResult}
\docType{data}
\title{
Rongelap island data 
}
\description{
A \code{SpatialPointsDataFrame} containing the Rongelap data, in a UTM projection.
}
\usage{
data("rongelapUTM")
data("rongelapResult")
}
\details{
These coordinates were obtained by translating and rotating the original Rongelap data
until all the coordinates fit into the Rongelap border given by www.gadm.org.  
So they are not exact.
}
\source{
See the help file for \code{\link[geoRglm]{rongelap}}, or 
\url{http://www.leg.ufpr.br/doku.php/pessoais:paulojus:mbgbook:datasets}
}
\examples{
data("rongelapUTM")
if(require("sp", quietly=TRUE)){
	plot(rongelapUTM)
}

\dontrun{
rongelapBorderLL = raster::getData("GADM", 
country="MHL",level=0) 
library("rgdal")
rongelapBorderUTM = spTransform(rongelapBorderLL, 
	CRS(proj4string(rongelapUTM)))
plot(rongelapBorderUTM, add=TRUE)
}

rongelapUTM$logOffset = log(rongelapUTM$time)
rongelapUTM$site = seq(1, length(rongelapUTM$time)) 
  
forBugs = glmmBUGS(
formula=count + logOffset ~ 1, family="poisson",
    data=rongelapUTM@data, effects="site", 
    spatial=rongelapUTM,
    priors=list(phisite="dgamma(100,1)")
    )

startingValues = forBugs$startingValues
startingValues$phi=list(site = 100)

source("getInits.R")

\dontrun{
rongelapResult = bugs(forBugs$ragged, getInits, 
  parameters.to.save = names(getInits()),
    model.file="model.bug", n.chain=2, n.iter=20, n.burnin=4, n.thin=2,
      program="winbugs", debug=TRUE,
      working.directory=getwd())
rongelapParams = restoreParams(rongelapResult, forBugs$ragged)
}

data("rongelapResult")
rongelapParams = restoreParams(rongelapResult)


checkChain(rongelapParams)

rongelapSummary=summaryChain(rongelapParams)
}



\keyword{datasets}
