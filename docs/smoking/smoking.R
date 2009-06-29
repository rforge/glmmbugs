# read in the popdata @
library(gdata)
popdata <- read.xls("test.xls", header=TRUE)

popdata$sex = factor(popdata$sex, levels=c(1,2), labels=c("M", "F"))

# read in DA 2001, in C:\Documents and Settings\luzhou\Paul
library(maptools)
adjdata01 = readShapePoly(fn="yyy2001")

library(spdep)
#adj01 = poly2nb(adjdata01, as.character(adjdata01$DA2001))

# just grab 1 region: EA_or_DA = 35230258
hamiltoncchs = popdata[popdata$CD_sub=="3525005",]
hamiltonSpatial = adjdata01[adjdata01$CSDUID == "3525005",]
library(spdep)
adj01 = poly2nb(hamiltonSpatial, as.character(hamiltonSpatial$DAUID))

#smpopdata = popdata[popdata$EA_or_DA == "35230258" ,]


library(glmmBUGS)
# change the response variable to the 0 and 1's:
hamiltoncchs$SMK_01a <- hamiltoncchs$SMK_01A - 1

## add in quadratic term agesq :

hamiltoncchs$agesq = hamiltoncchs$age^2

## for males only 

hamiltoncchsm = hamiltoncchs[hamiltoncchs$sex=="M",]

#forBugs = glmmBUGS(SMK_01a ~ age:sex, effects="EA_or_DA", 
#family="bernoulli", data=hamiltoncchs, spatial=adj01)


source("C:\\Documents and Settings\\luzhou\\My Documents\\glmmBUGS\\pkg\\glmmBUGS\\R\\addSpatial.R")
source("C:\\Documents and Settings\\luzhou\\My Documents\\glmmBUGS\\pkg\\glmmBUGS\\R\\getStartingValues.R")
source("C:\\Documents and Settings\\luzhou\\My Documents\\glmmBUGS\\pkg\\glmmBUGS\\R\\winBugsRaggedArray.R")
source("C:\\Documents and Settings\\luzhou\\My Documents\\glmmBUGS\\pkg\\glmmBUGS\\R\\writeBugsModel.R")

#data = hamiltoncchs
effects="EA_or_DA"
spatial = adj01
spatialEffect = "EA_or_DA"
family="bernoulli"

formula = SMK_01a ~ age:sex
# formula for quadratic terms: 

formula = SMK_01a ~ (age+agesq):sex


# for males in the model only: 
data =hamiltoncchsm
formula = SMK_01a ~ age+agesq
reparam =list(c(25, 625))
names(reparam)= "observations"
prefix=NULL

    data = getDesignMatrix(formula, data, effects)
    data = na.omit(data)
    covariates = attributes(data)$covariates
    observations = attributes(data)$response
    
     #if new col created, change it back for ragged
     
     
#        reparam = paste(prefix, reparam, sep="")
     
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

#popdata$ageCut = cut(popdata$age, c(0, seq(20, 60, by=5), Inf))
smkParams = restoreParams(smkResult, ragged$ragged)

    ontarioSummary = summaryChain(ontarioParams)












