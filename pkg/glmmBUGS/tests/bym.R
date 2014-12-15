havePackages = c(
    'diseasemapping'=require('diseasemapping', quietly=TRUE), 
    "spdep"=require('spdep', quietly=TRUE), 
    'R2OpenBUGS'=require('R2OpenBUGS', quietly=TRUE)
)

print(havePackages)

if(all(havePackages)){
  data('kentucky')
  larynxRates = structure(c(0, 0, 0, 0, 1e-06, 6e-06, 2.3e-05, 4.5e-05, 9.9e-05, 
          0.000163, 0.000243, 0.000299, 0.000343, 0.000308, 0.000291, 0.000217, 
          0, 0, 0, 1e-06, 1e-06, 3e-06, 8e-06, 1.3e-05, 2.3e-05, 3.5e-05, 
          5.8e-05, 6.8e-05, 7.5e-05, 5.5e-05, 4.1e-05, 3e-05), .Names = c("M_10", 
          "M_15", "M_20", "M_25", "M_30", "M_35", "M_40", "M_45", "M_50", 
          "M_55", "M_60", "M_65", "M_70", "M_75", "M_80", "M_85", "F_10", 
          "F_15", "F_20", "F_25", "F_30", "F_35", "F_40", "F_45", "F_50", 
          "F_55", "F_60", "F_65", "F_70", "F_75", "F_80", "F_85"))

kentucky = getSMR(kentucky, larynxRates, larynx,
    regionCode="County")

kAdjMat = spdep::poly2nb(kentucky,
    row.names=as.character(kentucky$County))

library('glmmBUGS')
forBugs = glmmBUGS(observed + logExpected ~ poverty,
    effects="County", family="poisson", 
    spatial=kAdjMat,
    modelFile='bym.txt',
    data=kentucky@data
)

startingValues = forBugs$startingValues

source("getInits.R")

# for some reason OpenBUGS won't run unless
  # all the starting values for RCountySpatial are zero
int2 = function() {
  res = getInits()
  res$RCountySpatial = rep(0, length(res$RCountySpatial))
  res
}

# find patrick's OpenBUGS executable file
if(Sys.info()['user'] =='patrick') {	 
  obExec = system(
      "find /store/patrick/ -name OpenBUGS",
      TRUE)
  obExec = obExec[length(obExec)]
} else {
  obExec = NULL
}


kResult = R2OpenBUGS::bugs(forBugs$ragged, 
  inits=int2,
  model.file = "bym.txt", 
  n.chain = 2,
  n.iter = 500, n.burnin = 10,
  parameters = names(int2()),
  n.thin = 10, OpenBUGS.pgm = obExec
)

kParams = restoreParams(kResult,
    forBugs$ragged)
kSummary = summaryChain(kParams)

kSummary$scalars[,c('mean', 'sd')]

pdf("checkChainBym.pdf",height=4,width=8)
checkChain(kParams, c("intercept", "SDCountySpatial"))
dev.off()

kentucky = mergeBugsData(kentucky, kSummary)

if(require('sp', quietly=TRUE)) {
  pdf("spplotBym.pdf",height=4,width=8)
  print(spplot(kentucky, "FittedRateCounty.mean"))
  dev.off()
}

}