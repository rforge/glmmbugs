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
#formula = SMK_01a ~ age+agesq
reparam =list(c(0.14, 0.32, 0.27, 0.13, 0.05))
names(reparam)= "observations"
prefix=NULL
mypriors = c("SDEA_or_DA" = "dunif(0, 25)", "SDEA_or_DASpatial"="dunif(0, 100)")


ageRange = range(hamiltoncchsm$age)

#internalKnots = c(6, 10, 14)  # 3 internal knots
internalKnots = c(32, 52)  # 3 internal knots

#Get basis for cubic splines, 4+3-1(intercept) = 6
library(splines)
thesplines = as.matrix(bs(hamiltoncchsm$age, knots=internalKnots, Boundary.knots=ageRange))
colnames(thesplines) = paste("s", colnames(thesplines), sep="")   
class(thesplines)="matrix"
hamiltoncchsm = cbind(hamiltoncchsm, thesplines) 

formula = SMK_01a ~ s1+ s2+ s3 + s4 + s5
#hamiltoncchsm$agebs = bs(hamiltoncchsm$age, knots = c(32,52))
#hamiltoncchsm$agesqbs = bs(hamiltoncchsm$agesq, knots = c(3000, 4500))


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
        prefix= attributes(ragged)$prefix, reparam =reparam, priors= mypriors)

    library(R2WinBUGS)     
    source("getInits.R")

## check whether this is unparam or reparam:     

       smkResult = bugs(ragged, getInits, 
       parameters.to.save = names(getInits()),
       model.file="model.bug", n.chain=3, n.iter=50, n.burnin=10, n.thin=2,
       program="WinBUGS", debug=T)     

#popdata$ageCut = cut(popdata$age, c(0, seq(20, 60, by=5), Inf))
smkParams = restoreParams(smkResult, ragged)

smkSummary = summaryChain(smkParams)
source("C:\\Documents and Settings\\luzhou\\My Documents\\newDiseaseMapping\\pkg\\diseasemapping\\R\\mergeBugsData.R")

smk = mergeBugsData(hamiltonSpatial, smkSummary, by.x= "DAUID")













