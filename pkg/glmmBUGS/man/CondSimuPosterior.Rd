\name{CondSimuPosterior}
\alias{CondSimuPosterior}

\title{
Simulate from the posterior of a geostatistical random effect on a grid
}
\description{
After having obtained a posterior sample of model parameters and spatial random effects at the data locations, 
this function will simulate from the posterior of the random effects on a grid of locations.
}
\usage{
CondSimuPosterior(params, locations.obs, xgrid = NULL, ygrid = NULL, gridSize = NULL,thin=1)
}

\arguments{
  \item{params}{
Posterior sample of parameters, from the function \code{restoreParams}
}
  \item{locations.obs}{
Locations of the observations, either the ragged array sent to BUGS or 
a matrix or SpatialPoints object containing the locations.
}
  \item{xgrid}{
Range of x values for the grid to be simulated on, or a vector of x values.  Defaults to the 
range of x values of the observed data.
}
  \item{ygrid}{
As xgrid for the y values.
}
  \item{gridSize}{
The size of the grid cells, must be specified unless xgrid and ygrid are provided as vectors. 
}
\item{thin}{use only one out of every 'thin' samples, defaults to 1 which uses every sample.}
}
\details{
Uses the function \code{CondSimu} in the \code{RandomFields} package.
}
\value{A list with entries
  \item{x}{vector of x coordinates of the grid cells}
  \item{y}{vector of y coordinates of the grid cells}
  \item{z}{Four dimensional array containing posterior samples, with dimensions being x, y, chain, and iteration.}

}
\author{
Patrick Brown
}

\seealso{
The Rongelap example in \code{\link{glmmBUGS}}
}
\examples{
}

