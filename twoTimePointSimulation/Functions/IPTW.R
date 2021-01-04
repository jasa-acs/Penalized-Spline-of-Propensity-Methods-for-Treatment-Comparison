###############################################################################
####dataInput--input dataset
###linear--linear or nonlinear in covariates
##propenVarListA0.Init--propensity model at first time point
####propenVarListA1.Init--propensity model at second time point
###treat.varnameA0--first treatment variable name
####treat.varnameA1--second treatment variable name
###outcome.varname--final outcome of interest
###original=1 for original dataset, 0 for bootstrap sample

IPTW = function(dataInput, linear, propenVarListA0.Init, propenVarListA1.Init,
                treat.varnameA0, treat.varnameA1, outcome.varname, original=1) { 
  tryCatch ( 
    {
      
      data=NULL
      
      if (original==1){
        
        data=dataInput
        
      } else {
        
        num00=which( (dataInput[,treat.varnameA0]==0) & (dataInput[,treat.varnameA1]==0) )
        num01=which( (dataInput[,treat.varnameA0]==0) & (dataInput[,treat.varnameA1]==1) )
        num10=which( (dataInput[,treat.varnameA0]==1) & (dataInput[,treat.varnameA1]==0) )
        num11=which( (dataInput[,treat.varnameA0]==1) & (dataInput[,treat.varnameA1]==1) )
        
        data=dataInput[c(sample(x=num00, size=length(num00),replace=T),
                         sample(x=num01, size=length(num01),replace=T),
                         sample(x=num10, size=length(num10),replace=T),
                         sample(x=num11, size=length(num11),replace=T)),]
        
        
      }
      
      treatIndA0 = data[, treat.varnameA0]  ###indicator of treated subjects
      treatIndA1 = data[, treat.varnameA1]  ###indicator of treated subjects
      Yobs = data[, outcome.varname]
      
      
      ######## MSM weight########################################################### 
      model1aF=formulaF(varList=c(1), y.name=treat.varnameA0)
      model1a=glm(model1aF, data=data, family="binomial") 
      temp1a=treatIndA0 * model1a$fitted.values + (1-treatIndA0) * (1-model1a$fitted.values)
      
      model1bF=formulaF(varList=c("A0"), y.name=treat.varnameA1)
      model1b=glm(model1bF, data=data, family="binomial") 
      temp1b=treatIndA1 * model1b$fitted.values + (1-treatIndA1) * (1-model1b$fitted.values)
      
      numer=temp1a * temp1b 
      
      #########
      model2aF=formulaF(varList=propenVarListA0.Init, y.name=treat.varnameA0)
      model2a=glm(model2aF, data=data, family="binomial") 
      temp2a=treatIndA0 * model2a$fitted.values + (1-treatIndA0) * (1-model2a$fitted.values)
      
      data$L1L0 = data$L1-data$L0
      data$L1bL0b = data$L1b-data$L0b
      model2bF = formulaF(varList=propenVarListA1.Init, y.name=treat.varnameA1)
      model2b = glm(model2bF, data=data, family="binomial") ### correctly specified treatment mechanism 
      temp2b = treatIndA1 * model2b$fitted.values + (1-treatIndA1) * (1-model2b$fitted.values)
      
      demo=temp2a * temp2b 
      weight=numer/demo 
      
      ####################################################################################
      ### E(u11)
      term11=sum(weight * Yobs * as.numeric( (treatIndA0==1) & (treatIndA1==1) ) )/sum(weight * as.numeric((treatIndA0==1) & (treatIndA1==1)))
      
      ####################################################################################
      ### E(u10)
      term10=sum(weight * Yobs * as.numeric( (treatIndA0==1) & (treatIndA1==0) ) )/sum(weight * as.numeric((treatIndA0==1) & (treatIndA1==0)))
      
      ####################################################################################
      ### E(u01)
      term01=sum(weight * Yobs * as.numeric( (treatIndA0==0) & (treatIndA1==1) ) )/sum(weight * as.numeric((treatIndA0==0) & (treatIndA1==1)))
      
      ####################################################################################
      ### E(u00)
      term00=sum(weight * Yobs * as.numeric( (treatIndA0==0) & (treatIndA1==0) ) )/sum(weight * as.numeric((treatIndA0==0) & (treatIndA1==0)))
      
      estFinal_iptw11=term11-term00  ###estimate
      estFinal_iptw10=term10-term00  ###estimate
      estFinal_iptw01=term01-term00   ###estimate
      
      return( c(estFinal_iptw11, estFinal_iptw10, estFinal_iptw01) ) ###
      
    }, error=function(e) return(NA) )
  
}
