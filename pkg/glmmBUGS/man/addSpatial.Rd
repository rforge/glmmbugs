\name{addSpatial}
\alias{addSpatial}
\title{Calculate adjacency values for WinBUGS}
\description{
Put an adjacency object in a ragged array
}
\usage{
addSpatial(map, raggedArray, effect = NULL)
}

\arguments{
  \item{map}{ a spatialPolygonsDataFrame object, or an nb object or a list of two vectors, adj and num}
  \item{raggedArray}{ the result from winBugsRaggedArray }
  \item{effect}{ a character vector listing the effect names }
}
\details{
Computes the values need by the \code{car.normal} distribution in WinBUGS.
}
\value{
 The ragged array is returned, with the following additional elements
  \item{num}{a vector of the number of neighbours of each region}
  \item{adj}{a vector containing the neighbours}
  \item{weights}{a vector of ones, the same length as \code{adj}}
  \item{NregionSpatial}{where 'region' is replaced by the name of the effect.
  The number of regions.}
}
\references{Also see the geoBUGS manual}
\author{ Patrick Brown}
\seealso{   }
\examples{

\dontrun{
library(diseasemapping)
library(sp)
library(spdep)
data(popdata)
data(casedata)
model <- getRates(casedata, popdata, ~age*sex, S=c("M","F"))
popdatanew <- getSMR(model, popdata, casedata)
# if region ID's are stored as factors, make sure to convert
# them to characters rather than the default of converting them
# to integers
ontarioNB = poly2nb(popdata,row.names=as.character(popdata[["CDSUID"]]) )

popRagged = glmmBUGS(cases + expected ~ 1, data=popdatanew@data, effects="CSDUID",
  family="poisson", spatial=ontarioNB)
ragged = popRagged$ragged

raggedWithSpatial = addSpatial(popdatanew, ragged, "CSDUID")

writeBugsModel("model.bug", "CSDUID", NULL, c("count", "expected"), "poisson", spatial="CSDUID")
startingValues = popRagged$startingValues
source("getInits.R")
library(R2WinBUGS)
popResult = bugs(raggedWithSpatial, getInits, parameters.to.save = names(getInits()), model.file="model.bug", n.chain=3, n.iter=1000, n.burnin=100, n.thin=10, program="winbugs")   
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ }
