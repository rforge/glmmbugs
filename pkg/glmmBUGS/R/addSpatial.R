addSpatial = function(map, raggedArray=NULL, effect=NULL, prefix=NULL) {

  
  if(is.null(effect)) {
  # assume the effect is at the lowest level
  # the lowest effect should have an Neffect and Seffect component
    effect = grep("^(N|S)", names(raggedArray), value=T)
    effect = substr(effect, 2,  10000)
    effect = table(effect)
    effect = names(effect)[effect==2]
    if(length(effect)!= 1) {
      warning("not sure which effect is the spatial one")
      effect = "region"
    }
  }

  # compute the adjancy matrix, if map is an sp object
  if(any(slotNames(map)=="polygons")) {
    map = poly2nb(map,row.names=map[[effect]])
  } 
  if(class(map)=="nb") {
    theregions = attributes(map)$region.id
 #   map = list(num = sapply(map, length),
 #     adj = unlist(map)  )
	map=nb2WB(map)	
    names(map$num) = theregions

  } else {
    if(!all(c("adj", "num") %in% names(map))) {
      warning("adj or num are missing in the map")
    }
    theregions = names(map$num)
  }
  
  # names of the random effects
  effectNames = names(raggedArray[[paste("S", prefix, effect, sep="")]])
  effectNames = effectNames[effectNames != "end"]
  if(!length(effectNames))
    effectNames = theregions


  # check to see each effect name is a region
  if(!all(effectNames %in% theregions))
    warning("some random effects don't appear to be regions")


Sspatial = seq(1, length(theregions))
names(Sspatial) = theregions
raggedArray[[paste("Sspatial", prefix, effect, sep="")]] = Sspatial[effectNames]


  raggedArray[[paste("adj", prefix, effect,  sep="")]] = map$adj
  raggedArray[[paste("num", prefix, effect, sep="")]] = map$num
  if(is.null(map$weights)) {
    raggedArray[[paste("weights", prefix, effect, sep="")]]= rep(1, length(map$adj))
  } else {
    raggedArray[[paste("weights", prefix, effect, sep="")]] = map$weights
  }
  raggedArray[[paste("N", prefix, effect, "Spatial", sep="")]] = length(map$num)

  raggedArray

}
