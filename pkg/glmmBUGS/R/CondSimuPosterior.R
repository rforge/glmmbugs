CondSimuPosterior = function(params, locations.obs, xgrid=NULL, ygrid=NULL, gridSize=NULL) {
	
	thePhi = grep("^phi", names(params), value=T)[1]
	theEffect = gsub("^phi", "",thePhi)
	theSD = paste("SD", theEffect, sep="")
	
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
	result = array(NA, c(length(xgrid), length(ygrid), Nchain, Niter))

	library(RandomFields)
	
	for(Dchain in 1:Nchain){
		for(Diter in 1:Niter){

		result[,,Dchain, Diter]=CondSimu("S", given=locations.obs, 
				data=params$Rsite[Diter,Dchain,], 
    		x=xgrid, y=ygrid, grid=TRUE, model="exponential", 
    		param=c(mean=0, variance=params[[theSD]][Diter,Dchain]^2, nugget=0, 
					scale=params[[thePhi]][Diter,Dchain]), pch="")
	
		}
	}

	result = list(x=xgrid, y=ygrid, z=result)
#	return(list(result, stuff))
	return(result)
	
}