`getStartingValues` <-
function(pql, ragged) {
# make a vector of starting values
# observations and effects are character strings
# covariates is a covariate list given to winBugsRaggedArray
# pql comes from glmmPQLstrings

# fixed effects  

startingValues = list(
  intercept = pql$coef$fixed["(Intercept)"]  
)

covariates = pql$covariates
effects = pql$effects
observations = pql$observations

if(is.list(covariates)) {
  for(Deffect in names(covariates)) {
    if(!is.null(covariates[[Deffect]]))
      startingValues[[paste("beta",Deffect,sep="")]] =
        pql$coef$fixed[covariates[[Deffect]] ]
        
  }
} else {
   startingValues$betas = pql$coef$fixed[names(pql$coef$fixed) != "(Intercept)"] 
}

# random effects
# data in the same order as the ragged array
  subdata = pql$data[,pql$effects]
  if(length(pql$effects)==1) {
      subdata = encodeString(as.character(subdata),max(nchar(subdata)))
      theorder = order(subdata)
  } else {
  # convert numerics to characters to make the order compatible with ragged array
    for(D in pql$effects)  {
      subdata[,D] = encodeString(as.character(subdata[,D]),max(nchar(subdata[,D])))
    }
    theorder = do.call(order, subdata)
  }
  thedata=pql$data[theorder,]

for(Deffect in seq(length(pql$effects), 1)) {
  theE = pql$effects[Deffect]

  # get one row of data for each different value of the effect
  theS = ragged[[paste("S", theE, sep="")]]
  theS = theS[-length(theS)]
  thedata = thedata[theS,]

  # get predicted values
  thepred = predict(pql, newdata=thedata, level=Deffect)
  
  # strip white space to make sure everything's compatible
  names(theS) = gsub(" ", "", names(theS))
  names(thepred) = gsub(" ", "", names(thepred))
  
  # put them in the same order as the ragged array
  startingValues[[paste("R", theE, sep="")]] = thepred[names(theS)]  
  # set covariates at this level to zero
  thedata[,covariates[[theE]] ] = 0
}
#
## covariate matrix
startingValues$vars = lapply(pql$modelStruct$reStruct, function(x) pql$sigma^2 * as.matrix(x)) 
#

# spatial
spatialEffect = grep("^N[[:alnum:]]+Spatial$" , names(ragged), value=T)
if(length(spatialEffect) ) {
spatialEffect = paste("R", gsub("^N", "", spatialEffect), sep="")
  spatialEffectIndep = gsub("Spatial$", "", spatialEffect)
  spatialEffectIndepVar = gsub("^R", "", spatialEffectIndep)
  spatialEffectVar = paste(spatialEffectIndepVar, "Spatial", sep="")

  
for(D in 1:length(spatialEffect)) {
# starting value for the proportion of spatial effect
spatialFactor = 0.5

  thenames = names(startingValues[[spatialEffectIndep]])
  startingValues[[spatialEffect[D] ]] =
      pql$coef$random[[spatialEffectIndepVar]][
        thenames, ,drop=T ]*spatialFactor

  
  startingValues[[spatialEffectIndep[D]]] =
    startingValues[[spatialEffectIndep[D]]] -  startingValues[[spatialEffect[D]]]

  # if there are some regions without data, add zeros as their starting values
  startingValues[[spatialEffect[D] ]] = c(
      startingValues[[spatialEffect[D] ]],
      rep(0, ragged[[paste("N", spatialEffectVar[D], sep="")]] -
        ragged[[paste("N", spatialEffectIndepVar[D], sep="")]] )
      )

  # starting values for variances
    startingValues$vars[[spatialEffectVar[D] ]] =
    sqrt(
      startingValues$vars[[spatialEffectIndepVar[D] ]]^2*spatialFactor
      )

    startingValues$vars[[spatialEffectIndepVar[D] ]] = sqrt(
        startingValues$vars[[spatialEffectIndepVar[D] ]]^2*(1-spatialFactor)
      )

}

}


return(startingValues)

}

