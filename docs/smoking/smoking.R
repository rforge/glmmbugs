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

# change directory to billy:
#cchsSpatial = readShapePoly(fn="da2001_HamiltonCSD")

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


#formula = SMK_01a ~ age:sex
# formula for quadratic terms: 

#formula = SMK_01a ~ (age+agesq):sex

prefix=NULL
mypriors = NULL

#c("SDEA_or_DA" = "dunif(0, 100)", "SDEA_or_DASpatial"="dunif(0, 100)", 
#"betaobservations[2]" = "dunif(-10, 10)")


# for males in the model only: 

#internalKnots = c(6, 10, 14)  # 3 internal knots
internalKnots = 25 #c(32, 52)  # 3 internal knots
splineDegree=2
ageRange = range(hamiltoncchsm$age)

#Get basis for cubic splines, 4+3-1(intercept) = 6
library(splines)

bsLast = function(x, df=NULL, knots=NULL, degree=3, Boundary.knots=range(x)) {
  fromBs = bs(x, df, knots, degree, Boundary.knots, intercept=T)
  fromBs[,-dim(fromBs)[2] ]
}


#formula = SMK_01a ~ age+agesq
reparam = list(c(as.vector(bsLast(50, knots=internalKnots, Boundary.knots=ageRange, degree=splineDegree))))

#reparam =list(c(0.53, 0.39, 0.036, 0, 0))
names(reparam)= "observations"

thesplines = as.matrix(bsLast(hamiltoncchsm$age, 
  knots=internalKnots, Boundary.knots=ageRange, 
    degree=splineDegree))

colnames(thesplines) = paste("s", colnames(thesplines), sep="")   
class(thesplines)="matrix"


hamiltoncchsm = cbind(hamiltoncchsm, thesplines) 
data =hamiltoncchsm

formula = as.formula(paste("SMK_01a ~ ", paste(colnames(thesplines), collapse="+")))

#formula = SMK_01a ~ s1+ s2+ s3 + s4 + s5
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
       model.file="model.bug", n.chain=3, n.iter=200000, n.burnin=10000, n.thin=100,
       program="WinBUGS", debug=T)     

#popdata$ageCut = cut(popdata$age, c(0, seq(20, 60, by=5), Inf))

source("C:\\Documents and Settings\\luzhou\\My Documents\\glmmBUGS\\pkg\\glmmBUGS\\R\\restoreParams.R")
smkParams = restoreParams(smkResult, ragged)
t(apply(smkParams$betas, 3, function(qq) c(mean=mean(exp(qq)), probGTone=mean(qq>0), quantile(exp(qq), c(0.025, 0.975)))))


smkSummary = summaryChain(smkParams)
checkChain(smkParams)
source("C:\\Documents and Settings\\luzhou\\My Documents\\newDiseaseMapping\\pkg\\diseasemapping\\R\\mergeBugsData.R")

smksp = mergeBugsData(hamiltonSpatial, smkSummary, by.x= "DAUID")

#for cchs: 
cchssp = mergeBugsData(cchsSpatial, smkSummary, by.x="DA2001")

# exceed probability plot: 
exceed =  apply(smkParams$REA_or_DA, 3, function(qq) mean(qq > log(1.1)))
pop = mergeBugsData(smksp, exceed, "DAUID", newcol="prob10")
exceed =  apply(smkParams$REA_or_DA, 3, function(qq) mean(qq > log(1.2)))
pop = mergeBugsData(smksp, exceed, "DAUID", newcol="prob20")
exceed =  apply(smkParams$REA_or_DA, 3, function(qq) mean(qq > log(1.5)))
pop = mergeBugsData(smksp, exceed, "DAUID", newcol="prob50")

spplot(pop, "prob10", lty=0)

# do the fitted plot using the color control according to the slides:

Ncol = 7
library(classInt)
ci=classIntervals(smksp$FittedRateEA_or_DA.mean, Ncol, style="kmeans")
library(RColorBrewer)
colours = brewer.pal(Ncol, "YlOrRd")
colFac = findColours(ci, colours)
ontarioBox = list(xlim=c(-83, -73.72315), ylim=c(41.8, 48))
print(plot(smksp, col=colFac, lwd=0.2) )
legend("topright", fill=attr(colFac, "palette"), legend=names(attr(colFac, "table")), bg="white") 



# do matplot: 
#basisW<-matrix(poster[,,13:18], ncol=6, dimnames=list(NULL, c("s1W","s2W","s3W","s4W"
#,"s5W","s6W")))

b <- grep("^s", dimnames(smkParams$betas)[[3]], value=T)
bCoef= matrix(smkParams$betas[,,b], ncol=length(b), dimnames=list(NULL, b))

ageSeq = seq(min(hamiltoncchsm$age), max(hamiltoncchsm$age), by=0.1)

splinesToPlot = as.matrix(bsLast(ageSeq, knots=internalKnots, 
  Boundary.knots=ageRange, degree=splineDegree))

smkpredict = splinesToPlot %*% t(bCoef)
smkMean = apply(smkpredict, 1, mean)
smkQuantile = apply(smkpredict, 1, function(qq) quantile(qq, c(0.025, 0.975)))




matplot(ageSeq, smkQuantile[1,],  xlim =c(12, 95), ylim = c(-0.8, 6.5), xlab= "age", ylab= "smk", type="l", lty=1, col="green", main="smoking against age")
matlines(ageSeq, smkQuantile[2,], lty=1, col="green")
matlines(ageSeq, smkMean, lty=1)




