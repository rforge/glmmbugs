\name{addSpatial}
\alias{addSpatial}
\title{Calculate adjacency values for WinBUGS}
\description{
Put an adjacency object in a ragged array
}
\usage{
addSpatial(map, raggedArray, effect = NULL,prefix=NULL)
}

\arguments{
  \item{map}{ a spatialPolygonsDataFrame object, or an nb object or a list of two vectors, adj and num}
  \item{raggedArray}{ the result from winBugsRaggedArray }
  \item{effect}{ a character vector listing the effect names }
  \item{prefix}{Character string to be appended to variable names}
}
\details{
Computes the values need by the \code{car.normal} distribution in WinBUGS.  
This function is called by \code{\link{glmmBUGS}} when a spatial argument is provided, \code{addSpatial} is 
usually not called by a user.
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

# get a winbugs model and data ready, without a spatial effect
data(ontario)

forBugs = glmmBUGS(formula=observed + logExpected ~ 1,
  effects="CSDUID",   family="poisson",
  data=ontario)

# now add a spatial effect.  
# first, compute the adjacency matrix 
# if region ID's are stored as factors, make sure to convert
# them to characters rather than the default of converting them
# to integers
\dontrun{
#library(spdep)
#library(diseasemapping)
#data(popdata)
popDataAdjMat = poly2nb(popdata,row.names=as.character(popdata[["CSDUID"]]) )
}
data(popDataAdjMat)

# add the adjacency matrix to the ragged array  
raggedWithSpatial = addSpatial(popDataAdjMat, forBugs$ragged, "CSDUID")

# write a new bugs model with a spatial effect
writeBugsModel("model.bug", "CSDUID", NULL, c("count", "expected"), "poisson", spatial="CSDUID")
startingValues = forBugs$startingValues
source("getInits.R")
\dontrun{
library(R2WinBUGS)
popResult = bugs(raggedWithSpatial, getInits, parameters.to.save = names(getInits()), model.file="model.bug", n.chain=3, n.iter=1000, n.burnin=100, n.thin=10, program="winbugs")   
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ }
