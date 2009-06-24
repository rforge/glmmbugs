# read in the popdata @
library(gdata)
popdata <- read.xls("test.xls", header=TRUE)

popdata$sex = factor(popdata$sex, levels=c(1,2), labels=c("M", "F"))

# read in DA 2001, in C:\Documents and Settings\luzhou\My Documents\CRC
library(maptools)
adjdata01 = readShapePoly(fn="CRC_M9903_DA2001_QAIPPEadj")

library(spdep)
#adj01 = poly2nb(adjdata01, as.character(adjdata01$DA2001))

# just grab 1 region: EA_or_DA = 35230258
hamiltoncchs = popdata[popdata$CD_sub=="3525005",]
hamiltonSpatial = adjdata01[adjdata01$CSDUID == "3525005",]
adj01 = poly2nb(hamiltonSpatial, as.character(hamiltonSpatial$DAUID))

#smpopdata = popdata[popdata$EA_or_DA == "35230258" ,]


library(glmmBUGS)
hamiltoncchs$SMK_01a <- hamiltoncchs$SMK_01A - 1
forBugs = glmmBUGS(SMK_01a ~ age:sex, effects="EA_or_DA", 
family="bernoulli", data=hamiltoncchs, spatial=adj01)















#popdata$ageCut = cut(popdata$age, c(0, seq(20, 60, by=5), Inf))













