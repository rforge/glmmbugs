# read in the popdata @
library(gdata)
popdata <- read.xls("test.xls", header=TRUE)

popdata$sex = factor(popdata$sex, levels=c(1,2), labels=c("M", "F"))

# read in DA 2001, in C:\Documents and Settings\luzhou\My Documents\CRC
#library(maptools)
#adjdata01 = readShapePoly(fn="CRC_M9903_DA2001_QAIPPEadj")

#library(spdep)
#adj01 = poly2nb(adjdata01, as.character(adjdata01$DA2001))

# just grab 1 region: EA_or_DA = 35230258
hamiltoncchs = popdata[popdata$CD_sub=="3525005",]
#hamiltonSpatial = adjdata01[adjdata01$CSDUID == "3525005",]
#adj01 = poly2nb(hamiltonSpatial, as.character(hamiltonSpatial$DAUID))
# change the response variable to the 0 and 1's:
hamiltoncchs$SMK_01a <- hamiltoncchs$SMK_01A - 1


data = hamiltoncchs
formula = SMK_01a ~ age:sex
effects="EA_or_DA"
spatial = adj01
spatialEffect = "EA_or_DA"
family="bernoulli"
prefix = NULL
reparam = NULL
library(glmmBUGS)
source("C:\\Documents and Settings\\luzhou\\My Documents\\glmmBUGS\\pkg\\glmmBUGS\\R\\addSpatial.R")
source("C:\\Documents and Settings\\luzhou\\My Documents\\glmmBUGS\\pkg\\glmmBUGS\\R\\getStartingValues.R")
source("C:\\Documents and Settings\\luzhou\\My Documents\\glmmBUGS\\pkg\\glmmBUGS\\R\\winBugsRaggedArray.R")
source("C:\\Documents and Settings\\luzhou\\My Documents\\glmmBUGS\\pkg\\glmmBUGS\\R\\writeBugsModel.R")
source("C:\\Documents and Settings\\luzhou\\My Documents\\glmmBUGS\\pkg\\glmmBUGS\\R\\getDesignMatrix.R")   

    data = getDesignMatrix(formula, data, effects)
    data = na.omit(data)
    covariates = attributes(data)$covariates
    observations = attributes(data)$response
    
   
     
    ragged = winBugsRaggedArray(data, effects = effects, covariates = covariates, 
        observations = observations, prefix= prefix , reparam=reparam)  
    #ragged$Xreparam = ????
    if (!is.null(spatial)) 
        ragged = addSpatial(spatial, ragged, spatialEffect, prefix= prefix)
    thepql = glmmPQLstrings(effects = effects, covariates = covariates, 
        observations = observations, data = data, family = family)
    startingValues = getStartingValues(pql = thepql, ragged = ragged, prefix=prefix, reparam = reparam )
    
    startingFunction(startingValues, file="getInits.R")

    spatialEffect = grep("^N[[:graph:]]+Spatial$", names(ragged), 
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
    

       smkResult = bugs(ragged, getInits, 
       parameters.to.save = names(getInits()),
       model.file="model.bug", n.chain=3, n.iter=50, n.burnin=10, n.thin=2,
       program="WinBUGS", debug=T)