CondSimuPosterior = function(params, locations.obs, xgrid=NULL, ygrid=NULL, gridSize=NULL, thin=1) {
	
	
	if(any(slotNames(locations.obs)=="proj4string")) {
		theproj4string = locations.obs@proj4string
	} else {
		theproj4string=NULL
	}
	
	thePhi = grep("^phi", names(params), value=T)[1]
	theEffect = gsub("^phi", "",thePhi)
	theSD = paste("SD", theEffect, sep="")
	theEffectR = paste("R", theEffect, sep="")
	
	
	if(all(paste(c("x", "y"), "Spatial", theEffect, sep="")%in% names(locations.obs))) {
		locations.obs = cbind(x=locations.obs[[paste("xSpatial", theEffect, sep="")]],
				y=locations.obs[[paste("ySpatial", theEffect, sep="")]])
	}
	if(length(grep("^SpatialPoints", class(locations.obs))))
		locations.obs = locations.obs@coords

	
	
	if(is.null(xgrid))
		xgrid = range(locations.obs[,1])
	
	if(is.null(ygrid))
		ygrid = range(locations.obs[,2])
	
	if(length(xgrid)==2 & !is.null(gridSize)){
		if(is.null(gridSize))
			warning("specify gridSize or make xgrid and ygrid vectors")
		xgrid = seq(xgrid[1], xgrid[2], by=gridSize)
		
	}
	if(length(ygrid)==2 & !is.null(gridSize)){
		if(is.null(gridSize))
			warning("specify gridSize or make xgrid and ygrid vectors")
		ygrid = seq(ygrid[1], ygrid[2], by=gridSize)
		
	}
	
	Nchain = dim(params[[thePhi]])[2]
	Niter = dim(params[[thePhi]])[1]
	Siter = seq(from=1, to=Niter, by=thin)
	Siter2 = seq(1, length(Siter))
	names(Siter) = as.character(1:length(Siter))
	result = array(NA, c(length(xgrid), length(ygrid), Nchain, length(Siter)))

	library(RandomFields)
	
	for(Dchain in 1:Nchain){
		for(Diter in Siter2){

		result[,,Dchain, Diter]=
			CondSimu("S", given=locations.obs, 
				data=params[[theEffectR]][Siter[Diter],Dchain,], 
    		x=xgrid, y=ygrid, grid=TRUE, model="exponential", 
    		param=c(mean=0, variance=params[[theSD]][Siter[Diter],Dchain]^2, nugget=0, 
					scale=params[[thePhi]][Siter[Diter],Dchain]), pch="")
	
		}
	}

	result = list(x=xgrid, y=ygrid, z=result, proj4string=theproj4string)
	
#	return(list(result, stuff))
	return(result)
	
}