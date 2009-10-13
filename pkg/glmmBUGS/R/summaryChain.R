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
                 
   for(D in thenames[!themat]) {
    result[[D]] = t(apply(chain[[D]], 3, getRes))
   
   }
   
   return(result)

}
                                            