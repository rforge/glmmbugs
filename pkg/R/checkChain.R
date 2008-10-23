checkChain = function(chain, parameters=NULL) {
if(is.array(chain))
  chain = list(beta=chain)
if(is.null(parameters)) {
  betas = grep("^beta", names(chain), value=TRUE)
  betas = c( grep("intercept", names(chain), value=TRUE), betas)
  betas  = c(grep("^SD", names(chain), value=TRUE), betas)
 
  if(length(betas) > 100)
    betas = betas[1:100]
 


  # find out the number of parameters
  themat = unlist(lapply(chain[betas], is.matrix))

  thepars = c(betas[themat], 
    unlist(lapply(chain[betas[!themat]], function(qq)
      dimnames(qq)[[3]])))

  } else {
   thepars = parameters 
   betas="beta"
  }


  par(mfrow=c(ceiling(length(thepars)/4),min(c(4, length(thepars)))), mar=c(2,2,1,0))
  
  
  for(Dbeta in betas) {
     if(length(dim(chain[[Dbeta]])) == 2) {
        plotOne(chain[[Dbeta]], Dbeta)
     } else {
        thenames =dimnames(chain[[Dbeta]])[[3]] 
        for(Dpar in thenames[thenames %in% thepars] ) {
          plotOne(chain[[Dbeta]][,,Dpar], Dpar)
        }
     }   
  }


  
  invisible() 
}

plotOne = function(mat, main=NULL) {

 matplot(mat, lty=1, type="l", main=main)
}