
###final.mod-prediction model
###baselineVar-baseline covariates
###numCarlo-number of monte carlo runs use in G computation
###here I use A0-to denote the first treatment indicator
###treat-treatment indicator (for the treatment that you are interested in estimating the associated potential outcomes)

###treat.varname-variable name denoting first treatment
###in the simulation studies, A0-denotes treatment at first time point

gcomputeFunc=function(final.mod, data, baselineVar, numCarlo, treat.varname, outcome.varname, treat){
  
  ############################################################### 
  ### generate baseline covariates from empirical distributions########
  baselineSim=NULL
  for(ind in 1:length(baselineVar)){
    draw_g=NULL
    draw_g=cbind(draw_g, rnorm(numCarlo, mean=mean(data[,baselineVar[ind]]), sd=sd(data[,baselineVar[ind]])) )
    colnames(draw_g)=baselineVar[ind]
    baselineSim=cbind(baselineSim, draw_g)
    
  }
  baselineSim=data.frame(baselineSim, A0_g=rep(treat, numCarlo))
  names(baselineSim)[which(names(baselineSim) == "A0_g")]=treat.varname
  
  ####squared terms ###############
  baselineSim[ , "L1_sq"] = baselineSim[, "L1"]^2
  baselineSim[ , "L2_sq"] = baselineSim[, "L2"]^2
  baselineSim[ , "L3_sq"] = baselineSim[, "L3"]^2
  baselineSim[, "L1L2"] = baselineSim[, "L1"] * baselineSim[, "L2"]
  
  ###### generate the final outcome of interest Y 
  includeVar = NULL
  includeVar = names(final.mod$coef)[-c(1)]
  includeDesign=NULL
  for(ind in 1:length(includeVar)){
    includeDesign=cbind(includeDesign, baselineSim[, includeVar[ind]])
  }
  
  
  Y_mean=(cbind(rep(1,numCarlo), includeDesign)%*% final.mod$coef) 
  Y_g=rnorm(numCarlo, mean=Y_mean, sd=summary(final.mod)$sigma) 
  
  baselineSim = data.frame(baselineSim, Y=Y_g)
  names(baselineSim)[which(names(baselineSim) == "Y")] = outcome.varname
  return( baselineSim)
  
  
}




