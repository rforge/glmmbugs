library(glmmBUGS)
library(diseasemapping)
data(popdata)
data(casedata)
model = getRates(casedata, popdata, ~age*sex, breaks=seq(0, 90, by=10) )
ontario = getSMR(popdata,model,  casedata)
data = ontario@data
data$a = rpois(583, 10)
formula = observed + logExpected ~ a
#data(ontario)
#data=  ontario@data
effects="CSDUID"
data(popDataAdjMat)
spatial = popDataAdjMat
spatialEffect = "CSDUID"
family="poisson"
prefix="cancer"
### 
reparam = "CSDUID"
###   when we have individual level covariates:
#reparam=c("CSDUID","observations")



    if(length(prefix)){ 
    oldeffects = effects
    effects = paste(prefix, effects, sep="")
    }
    #if prefix is not NULL or not empty ,,create a new column
    if(!is.null(prefix) & length(strsplit(toString(prefix)," ")[[1]])){
       data[,effects] = data[,oldeffects]
        create=TRUE
    }

    if (!is.character(family)) {
        warning("family must be a character string, ie \"poisson\" or \"binomial\" ")
    }
    

    data = getDesignMatrix(formula, data, effects)
    data = na.omit(data)
    covariates = attributes(data)$covariates
    observations = attributes(data)$response
    
     #if new col created, change it back for ragged
     if(create){
       data[,oldeffects] = data[,effects]
       data[,effects]<-NULL
       effects = oldeffects
     }
     
        reparam = paste(prefix, reparam, sep="")
     
    ragged = winBugsRaggedArray(data, effects = effects, covariates = covariates, 
        observations = observations, prefix= prefix , reparam=reparam)  
    #ragged$Xreparam = ????
    if (!is.null(spatial)) 
        ragged = addSpatial(spatial, ragged, spatialEffect, prefix= prefix)
    thepql = glmmPQLstrings(effects = effects, covariates = covariates, 
        observations = observations, data = data, family = family)
    startingValues = getStartingValues(pql = thepql, ragged = ragged, prefix=prefix, reparam = reparam )
    
    startingFunction(startingValues, file="getInits.R")

    spatialEffect = grep("^N[[:alnum:]]+Spatial$", names(ragged), 
        value = T)
    spatialEffect = gsub("^N", "", spatialEffect)
    spatialEffect = gsub("Spatial$", "", spatialEffect)
    
    effects = paste(prefix, effects, sep="")

    
    writeBugsModel(file="model.bug", effects = effects, covariates = covariates, 
        observations = observations, family = family, spatial = spatialEffect,
        prefix= attributes(ragged)$prefix, reparam =reparam)

    library(R2WinBUGS)     
    source("getInits.R")

## check whether this is unparam or reparam:     
     if(length(reparam)) reparamIntercept<-paste("interceptUnparam",prefix,sep="")

       ontarioResult = bugs(ragged, getInits, 
       parameters.to.save = c(names(getInits()), reparamIntercept),
       model.file="model.bug", n.chain=3, n.iter=50, n.burnin=10, n.thin=2,
       program="WinBUGS", debug=T)




    ontarioParams = restoreParams(ontarioResult, ragged$ragged)

    ontarioSummary = summaryChain(ontarioParams)



































