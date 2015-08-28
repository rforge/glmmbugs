\name{winBugsRaggedArray}
\alias{winBugsRaggedArray}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Ragged Arrays for multilevel models in BUGS}
\description{
  Suitable for unbalanced data.
}
\usage{
winBugsRaggedArray(data, effects = names(data)[-length(names(data))], 
covariates = NULL, observations = names(data)[length(names(data))], 
returnData = FALSE,
prefix=NULL, reparam=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{A data frame containing the response, covariates, and group membership.}
  \item{effects}{A vector of character strings containing the grouping levels, from most general to most specific.  Defaults to the column names of data, excluding the last column.}
  \item{covariates}{A list with names corresponding to effects and each element being a vector of covariates applicable at that level}
  \item{observations}{A character string giving the column of observations, or a vector where the first element is the observations and the remaning are offsets.}
  \item{returnData}{If true, returns the re-ordered data frame as well as the data frame}
  \item{prefix}{Character string to be appended to variable names}
  \item{reparam}{Vector of effect names, reparametrize the intercept by subtracting the mean of covariates at this level.}
}
\details{
   This funciton creates a list of data suitable for passing to the \code{\link[R2WinBUGS:bugs]{bugs}} function, suitable for implimentation as a ragged array.  The output can be passed to \code{\link{getStartingValues}} to manipulate the output from \code{\link{glmmPQLstrings}}, and to \code{\link{restoreParams}} to restore the original parametrisation from bugs output.
   

}
\value{
A list with the following components
\item{Nxx}{The number of levels in the most general groupping}
\item{Syy}{Indexing sequences, one for each level. If yy is level n, level n+1 has elements Syy[1] to Syy[2]-1 belonging to the first category of level n.}
\item{Xyy}{Matrix or vector of covariates belonging to level yy}
\item{}{vector of observations.}

}
\references{"Handling unbalanced datasets" in the 
"Tricks: Advanced Use of the BUGS Language" section of the bugs manual, at 
\url{http://www.openbugs.net/Manuals/Tricks.html}}
\author{Patrick Brown, \email{patrick.brown@utoronto.ca}}
 
\seealso{\code{\link[R2WinBUGS:bugs]{bugs}} }
\examples{
library(nlme)
data(Muscle)
muscleRagged = winBugsRaggedArray(Muscle, effects="Strip",  
observations="conc",
  covariates=list(observations="length"))

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.

