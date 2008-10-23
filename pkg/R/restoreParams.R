`restoreParams` <-
function(bugsResult, ragged=NULL) {

thearray = bugsResult$sims.array
parnames = dimnames(thearray)[[3]]
# vector valued parameters
vecPars = grep("\\[[[:digit:]]+\\]$", parnames, value=TRUE)
# matrix valued parameters
matPars =  grep("[[:digit:]+],[[:digit:]]+\\]$", parnames, value=TRUE)
# scalar parameters
scPars = parnames[! parnames %in% c(vecPars, matPars)]


result = list()


# if there are any random slopes (matrix parameters), put the slopes in their
# own element of the result list
if(length(matPars)) {
# find the precision matrix
precisionIndex = grep("^T", matPars)
precisions = matPars[precisionIndex]
matPars = matPars[-precisionIndex]
result$precision = thearray[,,precisions]

# convert precisions to variances
precisions = unique(gsub("\\[[[:digit:]]+,[[:digit:]]+\\]", "", precisions))


for(D in precisions) {
  result[[paste("var", D, sep="")]] = cholInvArray(result$precision, D)
}

# find the last column, which is the intercept
colno = substr(matPars, regexpr(",[[:digit:]]+\\]",  matPars)+1, 10000)
colno = substr(colno, 1, nchar(colno)-1)
maxcol = max(as.integer(colno))
if(is.na(maxcol))
  warning("can't find max col number")
# put the slopes in result
interceptcols = grep(paste(",", maxcol, "\\]", sep=""), matPars)
slopecols = matPars[-interceptcols]
result$slopes = thearray[,,slopecols]
# turn the intercept columns into vector parameters
matPars = which(parnames %in% matPars)
parnames[matPars] = gsub(paste(",", maxcol, "\\]", sep=""), "\\]", parnames[matPars])
dimnames(thearray)[[3]] = parnames
# update vecPars to reflect these intercepts
vecPars = grep("\\[[[:digit:]]+\\]$", parnames, value=TRUE)
} 


if(!length(scPars))
  warning("no parameter names")


for(D in scPars)
  result[[D]] = thearray[,,D]


# find grouping variables, all the variables with one dimensional indices
groups = unique(gsub("\\[[[:digit:]]+\\]$", "", vecPars))

# return(list(thearray, groups))

for(D in groups) {
  thisGroup = grep(paste("^", D, "\\[", sep=""), vecPars, value=TRUE)
  result[[D]] = thearray[,,thisGroup]
}

#return(list(result, groups))
    
  # if a ragged option is given, 
  # undo the reparametrisation and add better names to parameters
  if(!is.null(ragged)) {
     groups = paste("S", substr(groups, 2, nchar(groups)), sep="")

     randomEffects = groups[groups %in% names(ragged)]
     # order random effects 
     randomEffects = names(sort(unlist(lapply(ragged[randomEffects], length))))
     randomEffects = substr(randomEffects, 2, nchar(randomEffects))

 
     if(!length(randomEffects))
      warning(paste(toString(groups), ": can't find random effects in ragged object"))
      
     themeanOld = array(result[["intercept"]], 
        c(dim(result[["intercept"]]), 1))
     Nchain = dim(themeanOld)[2]
     torep = rep(1, length(ragged[[paste("S", randomEffects[1], sep="")]])-1)
     
#     return(list(result, groups, randomEffects))
     
     for(D in randomEffects) {  
        theS = ragged[[paste("S", D, sep="")]]
        theR = paste("R", D, sep="")
        thenames = names(theS)[-length(theS)]
   
        if(length(thenames) != (dim(result[[theR]])[3]) )
          warning(D, "different dimensions in bugsResult and ragged")
       dimnames(result[[theR]])[[3]] = thenames   
       
       DX =  paste("X", D, sep="")
       Dbeta = paste("beta", D, sep="")

      

       themean = themeanOld[,,torep]
       
       # expand the previous mean vector to the number of current random effects
#       return(list(themean=themean, rag = ragged[[DX]], res = result[[Dbeta]][,,]))
      if(!is.null(ragged[[DX]])) {
       # if there are covariates at this level
       theX = t(ragged[[DX]])
       thebeta = result[[Dbeta]]
       if(is.matrix(thebeta))
        thebeta = array(thebeta, c(dim(thebeta), 1))
       for(Dchain in 1:Nchain) {
          themean[,Dchain,] = themean[,Dchain,] + 
            thebeta[,Dchain,] %*% theX
       }  
      } 

      themeanOld = result[[theR]]
       result[[theR]] = result[[theR]] - themean
       
       
       torep = diff(theS)
       torep = rep(1:length(torep), torep) 

 
     }
     
     
     
     fixedEffects = grep("^X", names(ragged), value=TRUE)
     fixedEffects = substr(fixedEffects, 2, nchar(fixedEffects))
     
     Sbeta = paste("beta", fixedEffects, sep="")

     SX = paste("X", fixedEffects, sep="")

#     return(list(result, fixedEffects, SX, Sbeta))

     for(D in 1:length(fixedEffects)) {
        if(is.matrix(ragged[[SX[D]]])) {
          newnames = dimnames(ragged[[SX[D]]])[[2]]
          newnames = newnames[!(paste("beta", newnames, sep="") %in% names(result) )]
          if(length(newnames) == (dim(result[[Sbeta[D]]])[3]) )
            dimnames(result[[Sbeta[D]]])[[3]] = newnames
        }
     } 
     
     # make one matrix of betas
     betanameIndex = grep("beta", names(result))
      betas = NULL
     library(abind)
     for(D in betanameIndex) {
        if(is.matrix(result[[D]])) {
          result[[D]] = array(result[[D]], c(dim(result[[D]]), 1), 
            dimnames = list(dimnames(result[[D]])[[1]], dimnames(result[[D]])[[2]], 
            sub("beta", "", names(result)[D])) 
          )
        }
          betas = abind(betas, result[[D]])
     }
     result = result[-betanameIndex]
     result$betas = betas
    
  
  }
  
  
  return(result)
}

