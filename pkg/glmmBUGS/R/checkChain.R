checkChain = function(chain, parameters=NULL, oneFigure=TRUE) {

	if(is.array(chain))
  chain = list(beta=chain)

if(is.null(parameters)) {
	parameters = grep("^sd|range|intercept|betas",
			names(chain),value=TRUE,ignore.case=TRUE)
}

if(!all(parameters %in% names(chain)))
	warning("parameter ", parameters[!parameters %in% names(chain)], 
			"not found")

  # find out the number of parameters
  	scalars = unlist(lapply(chain[parameters], is.matrix))
	notScalars = parameters[!scalars]
	scalars = parameters[scalars]


	
	if(oneFigure) {
		Nplots = length(scalars)

		for(Dbeta in notScalars)
			Nplots = Nplots + dim(chain[[Dbeta]])[3]

	  par(mfrow=c(ceiling(Nplots/4),min(c(4, Nplots))))
  }
  
  for(Dbeta in parameters) {
	  if(Dbeta %in% scalars) {
		  plotOne(chain[[Dbeta]], Dbeta)
	  } else {
		  for(D2 in dimnames(chain[[Dbeta]])[[3]])
			  plotOne(chain[[Dbeta]][,,D2], 
					  D2)
		  }
	 }
  
  invisible() 
}

plotOne = function(mat, main=NULL) {
 matplot(mat, lty=1, type="l", xlab='iteration', ylab=main)
}