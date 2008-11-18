`writeBugsModel` <-
function(file, effects, covariates, observations, 
  family=c("bernoulli", "binomial", "poisson", "normal",  "other"),
  spatial=NULL) {

# spatial is a character string of names of random effects

  if(!is.character(family)) {
    warning("family must be a character string, ie \"poisson\" or \"binomial\" ")
  }
  family = family[1]
  if(!all(spatial %in% effects))
    warning("spatial effects used which are not specified as random effects")
  
  offset = observations[-1]
  observations = observations[1]
  
  if(family=="poisson") {
    link = "log"
    ddist = "dpois"
  } else if(family=="bernoulli") {
    link="logit"
    ddist="dbern"
  } else if(family=="binomial") {
    link="logit"
    ddist="dbin"
  } else if(family %in% c("normal", "gaussian")) {
    link=""
    ddist="dnorm"
  } else {
    ddist = paste("d", family, sep="")
    link=""
  }
  if(link != "") {
    link = paste(link, "(", sep="")
    endlink = ")"
  } else {
    endlink=""
  }
  
  

  sink(file)
  
  cat("model{\n\n")
  
  # the first effect
  Deffect=1
  indent = encodeString(" ", width=2*Deffect)
  theE = effects[Deffect]
  theD =  paste("D", theE, sep="") 
  cat("for(", theD, " in 1:N", theE, ") {\n\n", sep="")
  cat(indent, "R", effects[Deffect], "[", theD, "] ~ dnorm(mean", 
   theE, "[", theD, "], T", theE, ")\n",sep="")
  cat(indent, "mean", theE, "[", theD, "] <- intercept", sep="")
  
  # the covariates
  # check to see if there's more than one
  if(length(covariates[[theE]])==1) {
    cat(" + beta", theE, " * X", theE, "[", theD, "]", sep="")
  } else if (length(covariates[[theE]]) > 1) {
    cat(" + inprod2(beta", theE, "[] , X", theE, "[", theD, ",])", sep="")
    
  }   
  # spatial
  if(theE %in% spatial) {
     cat("+ R", theE, "Spatial[Sspatial", theE, "[", theD, "]]", sep="")
  }
  cat("\n")
  
  # subsequent effects, if any
  if(length(effects)>1) {
    for(Deffect in seq(2, length(effects))) {
      theE = effects[Deffect]
      theD = paste("D", theE, sep="")
      thePastD = paste("D", effects[Deffect-1], sep="")
      thePastS = paste("S", effects[Deffect-1], sep="")
      cat("\n", indent, "for(", theD, " in ", thePastS, "[",
        thePastD, "]:(", thePastS, "[", thePastD, "+1]-1)){\n\n", sep="")
      indent = encodeString(" ", width=2*Deffect)
      cat(indent, "R", theE, "[", theD, "] ~ dnorm(mean", 
        theE, "[", theD, "], T", theE, ")\n",sep="")
      cat(indent, "mean", theE, "[", theD, "] <- R", 
        effects[Deffect-1], "[",
        thePastD, "]", sep="")
  
  # the covariates
  # check to see if there's more than one
      if(length(covariates[[theE]])==1) {
        cat(" + beta", theE, " * X", theE, "[", theD, "]", sep="")
      } else if (length(covariates[[theE]]) > 1) {
        cat(" + inprod2(beta", theE, "[] , X", theE, "[", theD, ",])", sep="")
    
      }   
  # spatial
  if(theE %in% spatial) {
     cat("+ R", theE, "Spatial[", theD, "]", sep="")
  }

  cat("\n")
        
    }
    }

  # the observations
      theE = "observations"
      theD = paste("D", theE, sep="")
      thePastD = paste("D", effects[length(effects)], sep="")
      thePastS = paste("S", effects[length(effects)], sep="")
      
      # the loop
      cat("\n", indent, "for(", theD, " in ", thePastS, "[",
        thePastD, "]:(", thePastS, "[", thePastD, "+1]-1)){\n\n", sep="")
      indent = encodeString(" ", width=2*length(effects)+2)
      # distribution of observations
      cat(indent, observations, "[", theD, "] ~ ", ddist, "(mean", 
        theE, "[", theD, "]", sep="")
      # if binomial, add offsets
      if( family=="binomial" & length(offset)) {
        cat(", ")
        cat(toString(paste(offset, "[", theD, "]", sep="")))  
      } else if(family %in% c("normal", "gaussian"))
        cat(", Tobservations")  
      cat(")\n",sep="")
      # mean of observations  
      cat(indent, link, "mean", theE, "[", theD, "]", endlink, " <- R",  
        effects[length(effects)], "[", thePastD, "]", sep="")
    
      if(length(covariates[[theE]])==1) {
        cat(" + betaobservations * Xobservations[", theD, "]", sep="")
      } else if (length(covariates[["observations"]]) > 1) {
        cat(" + inprod2(betaobservations[] , Xobservations[", theD, ",])", sep="")
      }
      if(family!="binomial" & length(offset)) {
           cat("+", gsub(",", "+", toString(paste(offset, "[", theD, "]", sep="")))) 
      }   
  cat("\n\n")
    
    
    # the closing brackets
    effects = c(effects, "observations")
    
    for(Deffect in seq(length(effects),1)) {
       cat(encodeString("", width=2*Deffect-2), "}#", effects[Deffect],"\n",sep="")
    }
  
  # the spatial distributions
  for(Deffect in spatial) {  
    cat("R", Deffect, "Spatial[1:N", Deffect, "Spatial] ~ car.normal(adj", Deffect, 
        "[], weights", Deffect, "[], num", Deffect, "[], T", Deffect, "Spatial)\n", sep="")  
  }  
  
  # the priors
  cat("\n\n# priors\n\n")
  cat("intercept ~ dflat()\n")
  for(Deffect in effects) {
    thiscov = covariates[[Deffect]]
    if(length(thiscov)==1) {
        cat("beta", Deffect, " ~ dflat()\n", sep="")
    } else if(length(thiscov)>1) {
      for(Dpar in 1:length(covariates[[Deffect]]))
          cat("beta", Deffect, "[", Dpar,  "] ~ dflat()\n",sep="")
    }
  }
  cat("\n")
  if(! family %in% c("normal", "gaussian"))
    effects = effects[-length(effects)]
  for(Deffect in effects) {
    cat("T", Deffect, " <- pow(SD", Deffect, ", -2)\n", sep="")
    cat("SD", Deffect, " ~ dunif(0, 100)\n", sep="")
  }
  for(Deffect in spatial) {
       cat("T", Deffect, "Spatial <- pow(SD", Deffect, "Spatial, -2)\n", sep="")
       cat("SD", Deffect, "Spatial ~ dunif(0, 100)\n", sep="")
  }
  
  cat("\n} # model\n") 

  sink()
}

