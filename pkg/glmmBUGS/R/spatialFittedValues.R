spatialFittedValues = function(chain,threashold=1) {
	
	result=list()
	
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
				
	
				expWithInt = exp(chain[[D]]$z + array(chain$intercept, dim(chain[[D]]$z)) )
				
				 themean = apply(expWithInt, 1:2, mean)
				
				 thelist$z = themean 
					result[[D]]=list()
					
					result[[D]]$mean = raster(thelist, crs=theCRS)
					
					ssq = apply(expWithInt, 1:2, function(qq) sum(qq*qq))
					
					thesd = sqrt(ssq/prod( dim(chain[[D]]$z)[3:4]) - themean*themean)
					
					thelist$z = thesd 
					
					result[[D]]$sd = raster(thelist, crs=theCRS)
					
					thelist$z = apply(expWithInt>threashold, 1:2, mean)
					result[[D]]$probExc = raster(thelist, crs=theCRS)
				
				
			}
			
		} else{
			
			for(D in theGrids) {
				result[[D]] = list(
						x=chain[[D]]$x,y=chain[[D]]$y)
				expWithInt = exp(chain[[D]]$z + array(chain$intercept, dim(chain[[D]]$z)) )
				
				result[[D]]$mean= apply(expWithInt, 1:2, mean)
				
				result[[D]]$sd = apply(expWithInt, 1:2, sd)
				
				result[[D]]$probExc = apply(expWithInt>threashold, 1:2, mean)
			}
			
		} # end no raster	
	} # end if have grids   
	
	if(length(result)==1) result= result[[1]]

	result

}