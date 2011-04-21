\name{writeBugsModel}
\Rdversion{1.1}
\alias{writeBugsModel}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Write a bugs model file for a Generalised Linear Mixed Model
}
\description{
Given a list of effect groups, and the covariates associated with each level, a bugs model file is written using ragged arrays corresponding to output from \code{\link{winBugsRaggedArray}}        
}
\usage{
writeBugsModel(file, effects, covariates, observations, 
family = c("bernoulli", "binomial", "poisson", "normal", "other"), spatial = NULL, 
geostat=FALSE, 
prefix = "",reparam=NULL, brugs=FALSE, priors=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{file}{
a character string denoting the name of the bugs model file writen.
}
  \item{effects}{
vector of effect groups
}
  \item{covariates}{
A list with names corresponding to effects and each element being a vector of covariates applicable at that level
}
  \item{observations}{
A character string giving the column of observations, or a vector where the first element is the observations and the remaning are offsets.
}
  \item{family}{
Response distribution
}
  \item{spatial}{
name of the spatial random effect
}
\item{geostat}{Is this a geostatistical random effect?  Defaults to FALSE for the Besag, York and Mollie discrete spatial variation model}
  \item{prefix}{
the prefix 
}
  \item{reparam}{vector of random effect names, subtract covariates at this level from the intercept.
}
\item{brugs}{make the model file compatible with BRugs by using the inprod function in place of inprod2}
\item{priors}{character string of prior distributions, with the name of each element referring to the parameter it is the prior for}

}

\details{
The arguments to the function specify a generalised linear mixed model.  A file containing code for a corresponding bugs model is written.  The model uses ragged arrays to specify grouping factors, and includes covariates at the appropriate levels to aid in chain convergence.  It is intended that the user will edit this file before it's use.  The prior distributions in particular may not be appropriate.
}
\value{
A file, suitable for passing to the \code{\link[R2WinBUGS:bugs]{bugs}} function in \code{R2WinBUGS}.
}
\references{
"Handling unbalanced datasets" in the "Tricks: Advanced Use of the BUGS Language" 
section of the bugs manual, at \url{http://mathstat.helsinki.fi/openbugs/data/Docu/Tricks.html}
}
\author{
Patrick Brown, \email{patrick.brown@utoronto.ca}
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~
\section{Warning }{You are strongly encouraged to modify the model file prior to using it.}
\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
writeBugsModel("model.bug", effects="Strip",  observations="conc",
  covariates=list(observations="length"),
  family="normal", priors=c(intercept="dunif(-10,10)") ) 
}

