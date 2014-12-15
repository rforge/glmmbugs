havePackages = c(
    'R2OpenBUGS'= require('R2OpenBUGS', quietly=TRUE)
)

print(havePackages)

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

if(all(havePackages)) {
  
  # find patrick's OpenBUGS executable file
  if(Sys.info()['user'] =='patrick') {	 
    obExec = system(
      "find /store/patrick/ -name OpenBUGS",
    TRUE)
    obExec = obExec[length(obExec)]
  } else {
    obExec = NULL
  }
  
	bacResult = R2OpenBUGS::bugs(bacrag$ragged, inits=getInits,
			model.file = "bacteria.txt", n.chain = 3,
			n.iter = 600, n.burnin = 10,
			parameters = names(getInits()),
			n.thin = 4, OpenBUGS.pgm = obExec)
	
	bacParams = restoreParams(bacResult,
			bacrag$ragged)
	bacsummary = summaryChain(bacParams)
	
	bacsummary$betas[,c('mean', 'sd')]

  pdf("checkChainOpenbugs.pdf",height=4,width=8)
	checkChain(bacParams, c("intercept", "SDID"),oneFigure=TRUE)
  dev.off()
}