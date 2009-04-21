`glmmBUGS` <-
function(formula, data, effects, 
  modelFile = "model.bug", initFile = "getInits.R", 
  family=c("bernoulli", "binomial", "poisson", "gaussian"),
spatial=NULL, spatialEffect = NULL) {

  family = family[1]
  if(!is.character(family)) {
    warning("family must be a character string, ie \"poisson\" or \"binomial\" ")
  }




# get design matrix
data = getDesignMatrix(formula, data, effects)

covariates = attributes(data)$covariates
observations = attributes(data)$response


# get rid of rows where there are NA's in covariates
# but NA's in response are ok.
if(length(covariates)) {
  data = data[!apply(data[,unlist(covariates), drop=F], 1, 
    function(qq) any(is.na(qq))), ]
}


 # create the ragged array
ragged = winBugsRaggedArray(data, effects =effects,
  covariates = covariates, observations = observations)

#add spatial to the ragged array
if(!is.null(spatial))
	ragged = addSpatial(spatial, ragged, spatialEffect)

# run pql model

# get rid of missing values
dataForPQL = na.omit(data)

thepql = glmmPQLstrings(effects=effects, covariates = covariates, 
  observations = observations, data=dataForPQL, family=family)
 
# extract properly formatted starting values from the pql
startingValues = getStartingValues(pql=thepql, ragged=ragged) 

# write a function to generate starting values
startingFunction(startingValues, initFile)

# create the model
spatialEffect = grep("^N[[:alnum:]]+Spatial$", names(ragged), value=T)
spatialEffect = gsub("^N", "", spatialEffect)
spatialEffect = gsub("Spatial$", "", spatialEffect)

writeBugsModel(modelFile, effects =effects,
    covariates = covariates, observations = observations,
     family=family, spatial = spatialEffect ) # add spatial            
     
return(list(ragged=ragged, startingValues = startingValues, pql=thepql))  
}

