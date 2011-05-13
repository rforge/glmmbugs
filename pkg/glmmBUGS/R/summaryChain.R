summaryChain = function(chain, probs = c(0.005, 0.025, 0.05, 0.5)) {
   thenames = names(chain)
  thenames = thenames[! thenames %in% "deviance"]
   
   themat = unlist(lapply(chain[thenames], is.matrix))
   
   
   probs = unique(sort(c(probs, 1-probs)))
   
   getRes = function(qq) {
      thep = mean(qq>0)
      thep = min(thep, 1-thep)
   
      c(mean = mean(qq), 
        pval = thep,
        sd = sd(c(qq)),
        quantile(qq, probs = probs))
   }

    result = list(scalars = matrix(NA, sum(themat), length(getRes(1)),
      dimnames = list(thenames[themat], names(getRes(1)) ) )  )
    
    
   for(D in thenames[themat]) {
     result$scalars[D,] = getRes(chain[[D]])
   }
                
   vectorParams = thenames[!themat]
   vectorParams = vectorParams[!vectorParams %in% grep("Grid$", vectorParams, value=T)]
   
   for(D in vectorParams) {
    result[[D]] = t(apply(chain[[D]], 3, getRes))
   
   }
   
   # if this is a spatial model, get summaries of exponentials of fitted values
   
   if(length(grep("Spatial$",names(result)))) {
     for(D in grep("^Fitted", names(result), value=T)) {
        result[[paste("FittedRate", gsub("^FittedR", "", D), sep="")]]=
                t(apply(exp(chain[[D]]), 3, getRes))

     } 
   }
   
   # geostatistical variables on a grid
	theGrids = grep("Grid$", names(chain),value=T)
	if(length(theGrids)){
		haveRaster = library(raster, logical.return=T)
		
	if(haveRaster){
		
		for(D in theGrids) {
			if(!is.null(chain[[D]]$proj4string)) {
				theCRS = chain[[D]]$proj4string
			} else {
				theCRS=""
			}
	
			thelist = list(x=chain[[D]]$x,y=chain[[D]]$y)
		
			thelist$z = apply(chain[[D]]$z, 1:2, mean)
			result[[D]]$mean = raster(thelist, crs=theCRS)
	
			ssq = apply(chain[[D]]$z, 1:2, function(qq) sum(qq*qq))
			
	thelist$z = sqrt(ssq/prod( dim(chain[[D]]$z)[3:4])-thelist$z*thelist$z ) 
	result[[D]]$sd = raster(thelist, crs=theCRS)
	
	thelist$z = apply(chain[[D]]$z>0, 1:2, mean )
	result[[D]]$pgt0 = raster(thelist, crs=theCRS)
	

	
}
	
	} else{
		
	for(D in theGrids)
		result[[D]] = list(
				x=chain[[D]]$x,y=chain[[D]]$y,
				mean = apply(chain[[D]]$z, 1:2, mean),
				sd = apply(chain[[D]]$z, 1:2, sd),
				pgt0 = apply(chain[[D]]$z, 1:2, function(qq) mean(qq>0)))

	
	}
} # end if have grids   
	return(result)

}
                                            