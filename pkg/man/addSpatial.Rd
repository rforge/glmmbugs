\name{addSpatial}
\alias{addSpatial}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Calculate adjacency values for WinBUGS}
\description{
  Computes the values need by the \code{car.normal} distribution in WinBUGS.  If necessary, the regions are reordered to be in the same order as the ragged array, with any regions without observations put last.
}
\usage{
addSpatial(map, raggedArray, effect = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{map}{ ~~Describe \code{map} here~~ }
  \item{raggedArray}{ ~~Describe \code{raggedArray} here~~ }
  \item{effect}{ ~~Describe \code{effect} here~~ }
}
\details{
}
\value{
 The ragged array is returned, with the following additional elements
  \item{num}{a vector of the number of neighbours of each region}
  \item{adj}{a vector containing the neighbours}
  \item{weights}{a vector of ones, the same length as \code{adj}}
  \item{NregionSpatial}{where 'region' is replaced by the name of the effect.  The number of regions.}
}
\references{Also see the geoBUGS manual}
\author{ ~~who you are~~ }
\note{ ~~further notes~~ 

 ~Make other sections like Warning with \section{Warning }{....} ~
}
\seealso{ ~~objects to See Also as \code{\link{help}}, ~~~ }
\examples{

\dontrun{
library(diseasemapping)
library(sp)
library(spdep)
data(popdata)
data(casedata)
model <- getRates(casedata, popdata, ~age*sex, S=c("M","F"))
popdatanew <- getSMR(model, popdata, casedata)
popRagged = glmmBUGS(cases + expected ~ 1, data=popdatanew@data, effects="CSDUID", family="poisson")
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
\keyword{ }% __ONLY ONE__ keyword per line
