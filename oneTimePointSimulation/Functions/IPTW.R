###############################################################################
###dataInput-the dataset
###propenVarlist-the variables in the propensity model
###treat.varname-names of the treatment 
###outcome.varname=name of the outcome of interest
##original-1 for calculating treatment effect estimate for the original dataset, instead of bootstrap samples

IPTW=function(dataInput, propenVarList, treat.varname, outcome.varname, original=1) {
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
      
      
      ######## propensity score model ########################################################### 
      ##fit the propensity model
      ### for this paper, I use logistic regression models 
      propen.model=formulaF(varList=propenVarList, y.name=treat.varname)

      
      ######## propensity score model ########################################################### 
      model2a=glm(propen.model, data=data, family="binomial", control = list(maxit = 50))
      temp2a=NULL
      temp2a=(treatInd==1)*predict(model2a, newdata = data, type = "response") + (treatInd==0)*(1-predict(model2a, newdata = data, type="response"))  ###estimate propensity of observed treatment
      weight=1/temp2a 

      #####stablized
      term1a_All=sum(weight * Yobs * (treatInd==1))/sum(weight * (treatInd==1))
      term1b_All=sum(weight * Yobs * (treatInd==0))/sum(weight * (treatInd==0))
      
      ATE_IPTW=term1a_All-term1b_All

      return(  ATE_IPTW ) ###
      
    }, error=function(e) return(NA) )
  
}
