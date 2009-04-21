\name{glmmBUGS}
\alias{glmmBUGS}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{A function to run Generalised Linear Mixed Models in Bugs}
\description{
  Creates ragged arrays, writes a model file, and generates sensible starting estimates.
}
\usage{
glmmBUGS(formula, data, effects, modelFile = "model.bug", initFile = "getInits.R", family = c("bernoulli", "binomial", "poisson", "gaussian"), spatial=NULL, spatialEffect = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
   \item{formula}{A formula for the fixed effects portion of the model} 
   \item{data}{A data frame containing the response, covariates, and group membership}
   \item{effects}{A vector of character strings containing the grouping levels, from most general to most specific}
   \item{modelFile}{File for saving the bugs model}
   \item{initFile}{File for saving the function for generating initial values}
   \item{family}{distribution of responses}
   \item{spatial}{Spatial adj matrix, in forms of ploygon, nb object or a list of adj and num}
   \item{spatialEffect}{ spatial variable from data } 
}
\details{


Consider the following model, where \eqn{Y_{ijk}}{Y_ijk} is the number of absences from individual k from class j in school k.  

\deqn{Y_{ijk} \sim  Poisson(\mu_i)}{Y_ijk ~ Poisson(mu_i)}
\deqn{\log(\mu_i) = \delta age_{ijk} \beta + classSize_{ij} \alpha + schoolCategory_i \gamma + U_i + V_{ij}}{log(mu_i) = intercept + age_ijk beta + classSize_ij alpha + schoolCategory_i gamma + U_i + V_ij }
\deqn{U_i \sim N(0, \sigma^2)}{U_i ~ N(0, sigma^2)}
\deqn{V_{ij} \sim N(0, \nu^2)}{V_ij ~ N(0, nu^2)}

Here there are covariates which apply to each of the three levels, and random effects at the school and class level. If \code{data} is a data frame with one line per individual, the following would impliment this model: 

\code{glmmBUGS(data, effects=c("school","class"), covariates = list(school="schoolCategory", class="classSize", observations="age"), observations = "absences"), family="poisson")}

To aid in convergence, the bugs model is actually the following:

\deqn{\log(\mu_i) = age_{ijk} \beta + V_{ij}}{log(mu_i) = age_ijk beta +  V_ij }

\deqn{V_{ij} \sim N(U_i + classSize_{ij} \alpha , \nu^2)}{V_ij ~ N(U_i + classSize_ij alpha , nu^2)}

\deqn{U_i \sim N(\delta + schoolCategory_i \gamma, \sigma^2)}{U_i ~ N(intercept + schoolCategory_i gamma, sigma^2)}

and the funciton \code{\link{restoreParams}} subtracts the means from the random effects to restore the original set of equations.

\code{glmmBUGS} calls the following functions: 
\describe{
\item{\code{\link{getDesignMatrix}}}{to convert factors and interactions to indicator variables and find which covariates apply at which levels}
\item{\code{\link{winBugsRaggedArray}}}{ to prepare the ragged array}
\item{\code{\link{glmmPQLstrings}}}{ estimate starting values}
\item{\code{\link{writeBugsModel}}}{ to create a model file}
\item{\code{\link{getStartingValues}} }{to extract starting values from the glmmPQL result}
\item{\code{\link{startingFunction}} }{to write a function to generate random starting values}
}
Type \code{glmmBUGS} on the R command line to see the source code, it provides a good summary of the roles of the various functions in the \code{glmmBUGS} package.


}
\value{
Returns a list with the ragged array, from  \code{\link{winBugsRaggedArray}}, and the list of starting values from \code{\link{getStartingValues}}.  Writes a model file and an initial value function.  Note that the initial value function in \code{initFile} will look for an object called \code{startingValues}, which does not exist as this is part of a list.  Either create \code{startingValues <- result$startingValues} or edit \code{initFile}.
}

\references{"Handling unbalanced datasets" in the "Tricks: Advanced Use of the BUGS Language" section of the bugs manual, at \url{http://mathstat.helsinki.fi/openbugs/data/Docu/Tricks.html\#HandlingUnbalancedDatasets}}
\author{Patrick Brown, \email{patrick.brown@utoronto.ca}}

\note{glmmBUGS uses the \code{inprod2} function, which isn't implimented in OpenBugs, the model file will have to be modified for use with OpenBUGS.}
\section{Warning}{You are strongly encouraged to modify the model file and the intial value function file prior to using them.}   

\seealso{\code{\link{winBugsRaggedArray}}, \code{\link{glmmPQLstrings}} , \code{\link{writeBugsModel}}, \code{\link{getStartingValues}}, \code{\link{startingFunction}},\code{\link[R2WinBUGS:bugs]{bugs}}}

\examples{ 
library(nlme)
data(Muscle)

muscleRagged = glmmBUGS(conc ~ length, data=Muscle, effects="Strip", family="gaussian")
startingValues = muscleRagged$startingValues

\dontrun{
source("getInits.R")
require(BRugs)
require(R2WinBUGS)  
muscleResult = bugs(muscleRagged$ragged, getInits, parameters.to.save = names(getInits()),
                model.file="model.bug", n.chain=3, n.iter=1000, n.burnin=100, n.thin=10,
                             program="openbugs") 
}
 
data(muscleResult)

muscleParams = restoreParams(muscleResult, muscleRagged$ragged)  
summaryChain(muscleParams)
checkChain(muscleParams)



# a spatial example
\dontrun{
library(diseasemapping)

data(popdata)
data(casedata)

model = getRates(casedata, popdata, ~age*sex)
ontario = getSMR(popdata, model, casedata)
ontario = ontario@data[,c("CSDUID","observed","logExpected")]

library(spdep)
popDataAdjMat = poly2nb(popdata,row.names=as.character(popdata[["CSDUID"]]))
}

data(popDataAdjMat)
data(ontario)

forBugs = glmmBUGS(formula=observed + logExpected ~ 1,
  effects="CSDUID",   family="poisson", spatial=popDataAdjMat,
  data=ontario)

startingValues = forBugs$startingValues

source("getInits.R")

\dontrun{
ontarioResult = bugs(forBugs$ragged, getInits, 
  parameters.to.save = names(getInits()),
    model.file="model.bug", n.chain=3, n.iter=50, n.burnin=10, n.thin=2,
      program="winbugs", debug=T)
}

data(ontarioResult)

ontarioParams = restoreParams(ontarioResult, forBugs$ragged)

ontarioSummary = summaryChain(ontarioParams)

# posterior probability of having 10x excess risk
postProb = apply(ontarioParams$FittedRateCSDUID, 3, function(x) mean(x>10))
hist(postProb)

\dontrun{

ontario = mergeBugsData(popdata, ontarioSummary)

spplot(ontario, "FittedRateCSDUID.mean")


ontario = mergeBugsData(ontario, postProb, newcol="postProb", by.x="CSDUID")
spplot(ontario, "postProb")



}
}

