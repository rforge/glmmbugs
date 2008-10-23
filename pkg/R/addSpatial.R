addSpatial = function(map, raggedArray=NULL, effect=NULL) {
  
  if(is.null(effect)) {
  # assume the effect is at the lowest level
  # the lowest effect should have an Neffect and Seffect component
    effect = grep("^(N|S)", names(ragged), value=T)
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
    library(sp)
    map = poly2nb(map,row.names=map[[effect]])
  }
  theregions = attributes(map)$region.id

  # names of the random effects
  effectNames = names(raggedArray[[paste("S", effect, sep="")]])
  effectNames = effectNames[effectNames != "end"]
  if(!length(effectNames))
    effectNames = theregions


  # check to see each effect name is a region
  if(!all(effectNames %in% theregions))
    warning("some random effects don't appear to be regions")

  map = nb2mat(map, style="B", zero.policy=T)
  mode(map) = "logical"

  # reorder the adjancy list if necessary
    newOrder = c(effectNames, theregions[ ! theregions %in% effectNames] )
    dimnames(map)[[2]] = dimnames(map)[[1]]
    map = map[newOrder, newOrder]


  raggedArray$adj = unlist(apply(map, 2, which))
  raggedArray$num = apply(map, 2, sum)
  raggedArray$weights = rep(1, length(raggedArray$adj))
  raggedArray[[paste("N", effect, "Spatial", sep="")]] = length(raggedArray$num)

  raggedArray

}
