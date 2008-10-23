`glmmBUGS` <-
function(formula, data, effects, 
  modelFile = "model.bug", initFile = "getInits.R", 
  family=c("bernoulli", "binomial", "poisson", "gaussian")) {

  family = family[1]
  if(!is.character(family)) {
    warning("family must be a character string, ie \"poisson\" or \"binomial\" ")
  }

# get design matrix
data = getDesignMatrix(formula, data, effects)

covariates = attributes(data)$covariates
observations = attributes(data)$response

 # create the ragged array
ragged = winBugsRaggedArray(data, effects =effects,
  covariates = covariates, observations = observations)

# run pql model
thepql = glmmPQLstrings(effects=effects, covariates = covariates, 
  observations = observations, data=data, family=family)

 
# extract properly formatted starting values from the pql
startingValues = getStartingValues(pql=thepql, ragged=ragged) 

# write a function to generate starting values
startingFunction(startingValues, initFile)

# create the model
writeBugsModel(modelFile, effects =effects,
    covariates = covariates, observations = observations,
     family=family)            
     
return(list(ragged=ragged, startingValues = startingValues, pql=thepql))  
}

