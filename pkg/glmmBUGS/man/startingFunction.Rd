\name{startingFunction}
\alias{startingFunction}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Write a function to generate random MCMC starting values}
\description{The code for the resulting function is saved in a file, 
to be edited and sourced in before calling WinBUGS.
}
\usage{
startingFunction(startingValues, file = "getInits.R")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{startingValues}{list returned from \link{getStartingValues}}
  \item{file}{character string giving the name of the file to write to}
}
\details{
 Given a list containing initial estimates of parameters and random effects, 
  a text file is produced containing code for a function to generate random 
  starting values for use with the bugs() function.  It is intended that the file 
  produced be checked and edited prior to use. 
}
\value{
  A file, with the name given by the 'file' argument, is written.
}
\references{}
\author{Patrick Brown, \code{patrick.brown@utoronto.ca}}
\note{}
\section{Warning }{You are strongly encouraged to edit the file to ensure the result is sensible} 

\seealso{\link{getStartingValues}, \code{\link[R2WinBUGS:bugs]{bugs}}}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.
sval = list(intercept=0, beta = 1:2, Rperson = rep(0, 5), vars=list(person=1))
startingFunction(sval)

}

