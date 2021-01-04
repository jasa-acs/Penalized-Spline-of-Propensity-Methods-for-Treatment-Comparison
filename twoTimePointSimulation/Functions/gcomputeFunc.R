
###inter.mod1-prediction model for the first intermediate outcome L1
###inter.mod2-prediction model for the second intermediate outcome L1b
###final.mod-prediction model for the final outcome of interest
###baselineVar-baseline covariates
###numCarlo-number of monte carlo simulations
###here I use A0-to denote the first treatment indicator
###A1-to denote the second treatment indicator
###treatmentSequence denote the potential treatment regimes (00, 01, 10, 11)

###treat.varnameA0-variable name denoting first treatment, treat.varnameA1-variable name denoting second treatment
###in the simulation studies, A0-denotes treatment at first time point, A1-denotes treatment at second time point

gcomputeFunc=function(inter.mod1, inter.mod2, final.mod, data, baselineVar, numCarlo, treat.varnameA0, treat.varnameA1,
                      treatSeq){
  
  ############################################################### 
  ### generate L0 from empirical distribution of L0 obtained from the data
  baselineSim=NULL
  for(ind in 1:length(baselineVar)){
    draw_g=NULL
    draw_g=cbind(draw_g, rnorm(numCarlo, mean=mean(data[,baselineVar[ind]]), sd=sd(data[,baselineVar[ind]])) )
    colnames(draw_g)=baselineVar[ind]
    baselineSim=cbind(baselineSim, draw_g)
    
  }
  baselineSim=data.frame(baselineSim, A0_g=rep(treatSeq[1], numCarlo))
  names(baselineSim)[which(names(baselineSim) == "A0_g")]=treat.varnameA0
  
  ############# generate first intermediate outcome L1   
  L1_mean=predict(inter.mod1, newdata = baselineSim)
  baselineSim=data.frame(baselineSim, L1=rnorm(numCarlo, mean=L1_mean, sd=summary(inter.mod1)$sigma) )
  
  ########### generate second intermediate outcome L1b 
  L1b_mean=predict(inter.mod2, newdata = baselineSim)
  baselineSim=data.frame(baselineSim, L1b=rnorm(numCarlo, mean=L1b_mean, sd=summary(inter.mod2)$sigma) )
  
  baselineSim=data.frame(baselineSim, A1_g=rep(treatSeq[2], numCarlo))
  names(baselineSim)[which(names(baselineSim) == "A1_g")]=treat.varnameA1
  baselineSim$L1L1bInt = baselineSim$L1 * baselineSim$L1b
  
  ###### generate the final outcome of interest Y 
  Y_mean=predict(final.mod, newdata = baselineSim)
  Y_g=rnorm(numCarlo, mean=Y_mean, sd=summary(final.mod)$sigma) 
  
  baselineSim=data.frame(baselineSim, Y=Y_g)
  names(baselineSim)[which(names(baselineSim) == "Y")]=outcome.varname
  
  return(baselineSim[, c(treat.varnameA0, treat.varnameA1, outcome.varname)])
  
  
}




