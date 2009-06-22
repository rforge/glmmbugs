

regionData = census2001[census2001$region=="hamilton",]
popdata = popdata[popdata$region=="hamilton",]

popdata$ageCut = cut(popdata$age, c(0, seq(20, 60, by=5), Inf))
popdata$sex = factor(popdata$sex, levels=c(1,2), labels=c("M","F"))


census2001adj = adjmat(regionData)


forBugs = glmmBUGS(SMK_081a ~ age:sex, effects="DA2001", family="bernoulli", data=popdata, spatial=census2001adj)



