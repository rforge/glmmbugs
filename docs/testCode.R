
library(MASS)
data(bacteria)
bacteria$group = substr(bacteria$ID, 1,2)

forBugs = glmmBUGS(????, effects = c("group", "ID"), ?????)




# get design matrix
data = getDesignMatrix(formula, data, effects)

# get rid of rows where there are NA's in covariates
# but NA's in response are ok.
data = data[, !apply(data[,covariates], 1, function(qq) any(is.na(qq))) ]

covariates = attributes(data)$covariates
observations = attributes(data)$response

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



# a spatial example

# change a few region's observed to NA

