library(MASS)
data(bacteria)
head(bacteria)

bacterianew <- bacteria
bacterianew$yInt = as.integer(bacterianew$y =="y")
levels(bacterianew$trt) <- c("placebo","drug", "drugplus")
library(glmmBUGS)
bacrag <- glmmBUGS(formula = yInt ~trt + week, data = bacterianew,effects = "ID",
		modelFile = "model.bug",family = "bernoulli")


  modelFile = "model.txt"; initFile = "getInits.R"; 
family = c("bernoulli", "binomial", "poisson", "gaussian"); 
spatial = NULL; spatialEffect = NULL; reparam = NULL; prefix=NULL; priors=NULL;
brugs=.Platform$OS.type=="unix"
formula = yInt ~trt + week; data = bacterianew;effects = "ID";
modelFile = "model.bug";family = "bernoulli"
reparam="ID"
pql=thepql; 


bacrag <- glmmBUGS(formula = yInt ~trt + week, data = bacterianew,effects = "ID",reparam="ID",
		modelFile = "model.bug",family = "bernoulli")

startingValues = bacrag$startingValues 
source("getInits.R")  
library(R2jags)
bacResult = jags(bacrag$ragged, getInits,model.file = "model.bug",
    	n.chain = 3,n.iter = 2000, n.burnin = 100,
		parameters.to.save = names(getInits()),n.thin = 10)


