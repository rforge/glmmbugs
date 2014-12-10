library('MASS')
data('bacteria')
bacterianew <- bacteria
bacterianew$yInt = as.integer(bacterianew$y == "y")
levels(bacterianew$trt) <- c("placebo",
    "drug", "drugplus")

library('glmmBUGS')
bacrag <- glmmBUGS(formula = yInt ~ trt + week, 
    data = bacterianew,
    effects = "ID", modelFile = "bacteria.txt",
    family = "bernoulli",brugs=TRUE)

source("getInits.R")
startingValues = bacrag$startingValues

if( require("R2jags", quietly=TRUE)) {
  bacResult = jags(bacrag$ragged, 
      inits=getInits,
      model.file = "bacteria.txt", n.chain = 3,
      n.iter = 600, n.burnin = 10,
      parameters = names(getInits()),
      n.thin = 4)
  
  bacParams = restoreParams(bacResult$BUGSoutput,
      bacrag$ragged)
  bacsummary = summaryChain(bacParams)
  
  bacsummary$betas[,c('mean', 'sd')]
  
  checkChain(bacParams, c("intercept", "SDID"),oneFigure=TRUE)
}