#########this script computes treatment effect estimate from g computation
###dataInput-the dataset
###treat.varname-names of the treatment 
###outcome.varname-name of the outcome variable
##outcomeVarList0-variables included in the outcome model for Y0
##outcomeVarList1-varialbes included in the outcome model for Y1
##original-1 for calculating treatment effect estimate on the original dataset, instead of bootstrap samples (for calculating standard error)
###numCarlo-number of monte carlo runs in g computation and estimating the augmentation term in AIPTW
###baseList-variables at baseline, draw from baseline empirical distributions (in the simulations there are 3 baseline covariates for one time point treatment)

gcompute= function(dataInput, linear, outcomeVarList0, outcomeVarList1, treat.varname, outcome.varname, original=1, numCarlo=2000) { 
  
  tryCatch ( 
    {
      
      data=NULL
      if (original==1){
        
        data=dataInput
        
      } else {
        
        ###stratified bootstrap based on treatment groups
        treatID=which(dataInput[, treat.varname]==1)
        controlID=which(dataInput[, treat.varname]==0)
        data = dataInput[c( sample(treatID,replace=T), sample(controlID,replace=T) ),]  ##stratified bootstraps
        
      }
      
      treatInd = data[, treat.varname]  ###indicator of treated subjects
      Yobs = data[, outcome.varname]

      #### use g-computation to compute beta(Qn) 
      ###prediction model for Y1 and Y0
      modely0 = formulaF(varList=outcomeVarList0, y.name=outcome.varname)
      modely1 = formulaF(varList=outcomeVarList1, y.name=outcome.varname)
      
      mod3_1=lm(modely1, data=data[treatInd==1,]) 
      mod3_0=lm(modely0, data=data[treatInd==0,]) 
      
      ## numCarlo number of runs in G-computation 
      ############################################################### 
      ###simulate potential outcomes under each treatment A0=1
      Gcomp1=gcomputeFunc(final.mod=mod3_1, data=data, baselineVar=baselineVar, numCarlo=numCarlo, treat.varname=treat.varname, outcome.varname=outcome.varname, treat=1)

      
      ############################################################### 
      ###simulate potential outcomes under each treatment A0=0
      Gcomp0=gcomputeFunc(final.mod=mod3_0, data=data, baselineVar=baselineVar, numCarlo=numCarlo, treat.varname=treat.varname, outcome.varname=outcome.varname, treat=0)

      
      estFinal=mean(Gcomp1[, outcome.varname])-mean(Gcomp0[, outcome.varname])
      return(estFinal)

      
    }, error=function(e) return(NA) )
  
} 

