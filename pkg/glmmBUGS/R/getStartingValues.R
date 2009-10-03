`getStartingValues` <-
function(pql, ragged, prefix=NULL, reparam=NULL) {
# make a vector of starting values
# observations and effects are character strings
# covariates is a covariate list given to winBugsRaggedArray
# pql comes from glmmPQLstrings

# fixed effects  

startingValues = list()  
startingValues[[paste("intercept",prefix, sep="")]]=  pql$coef$fixed["(Intercept)"]  


#names(startingValues)[1]<-paste(prefix,names(startingValues)[1],sep="")

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
#  thedata =pql$data
  # dont use observation level covariates in predictions
  thedata[,covariates$observations]=0
  
  
for(Deffect in seq(length(pql$effects), 1)) {
  theE = pql$effects[Deffect]
  theE = paste(prefix, theE, sep="")
  
  # get one row of data for each different value of the effect
  theS = ragged[[paste("S", theE, sep="")]]
  theS = theS[-length(theS)]

  thedata = thedata[theS,]
  # get predicted values
  thepred = predict(pql, level=Deffect, newdata=thedata)#[theorder]
  
  # strip white space to make sure everything's compatible
  names(theS) = gsub(" ", "", names(theS))
  names(thepred) = gsub(" ", "", names(thepred))
  
  # put them in the same order as the ragged array
  theseStartingValues = rep(0, length(theS) )
  names(theseStartingValues) = names(theS)
  
  theseStartingValues[names(thepred)] = thepred
  
    
  startingValues[[paste("R", theE, sep="")]] = theseStartingValues

  # set covariates at this level to zero
  thedata[,covariates[[theE]] ] = 0
}
#
## covariate matrix
names(pql$modelStruct$reStruct)<-paste(prefix,names(pql$modelStruct$reStruct),sep="")
startingValues$vars = lapply(pql$modelStruct$reStruct, function(x) pql$sigma^2 * as.matrix(x)) 
#

# spatial
spatialEffect = grep("^N[[:graph:]]+Spatial$" , names(ragged), value=T)
if(length(spatialEffect) ) {
spatialEffect = paste("R", gsub("^N", "", spatialEffect), sep="")
  spatialEffectIndep = gsub("Spatial$", "", spatialEffect)
  spatialEffectIndepVar = gsub("^R", "", spatialEffectIndep)
  spatialEffectVar = paste(spatialEffectIndepVar, "Spatial", sep="")

  
for(D in 1:length(spatialEffect)) {
# starting value for the proportion of spatial effect
spatialFactor = 0.5

  # create a vector of zeros for starting values for the spatial component
  # getting names from the Sspatial index in the ragged array
  theStart = rep(0, ragged[[gsub("^R", "N", spatialEffect[D])]])
  spatialSeqName = gsub("^R", "Sspatial", spatialEffect[D])
  spatialSeqName = gsub("Spatial$", "", spatialSeqName)
  names(theStart) = names(ragged[[spatialSeqName]])

  # replace those zeros with spatialFactor times the random effect
  # from the pql model
  
  # find the names of spatial regions included in the pql model
  thenames = rownames(pql$coef$random[[spatialEffectIndepVar[D] ]])
  #names(startingValues[[spatialEffectIndep[D] ]])
  thenames = names(theStart)[names(theStart) %in% thenames]
  theStart[thenames] =  pql$coef$random[[spatialEffectIndepVar[D] ]][
        thenames, ,drop=TRUE]
  
  startingValues[[spatialEffect[D] ]] = theStart
  
  startingValues[[spatialEffectIndep[D]]][thenames] =
    startingValues[[spatialEffectIndep[D]]][thenames] -  
      theStart[thenames]
  
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


 if(is.character(reparam))   {
reparamName = reparam
reparam = list()
for(D in reparamName)
reparam[[D]] = NULL
 } 

 for(D in names(reparam)){
     if(D %in% names(covariates)){
       theName= paste("X", D, "reparam", sep="")  
       startingValues[[paste("intercept",prefix, sep="")]] = startingValues[[paste("intercept",prefix, sep="")]] +sum(pql$coef$fixed[covariates[[D]]] * ragged[[theName]])
      } 
 }  

 
return(startingValues)

}

                                                       