#dataOR=simdat; propenVarList=propenVarList; outcomeVarList0=outcomeVarList0; 
#outcomeVarList1=outcomeVarList1;
#treat.varname=treat.varname; outcome.varname=outcome.varname;
#original=1
### in this paper, we used the simplest basis functions-linear bases
###dataOR=original dataset
###fit models on bootstrap samples
###propenVarList-propensity model
##outcomeVarList0-model for Y0
###outcomeVarList1-model for Y1
###numKnots-number of knots truncated linear basis functions

pencomp=function(dataOR, propenVarList, outcomeVarList0, outcomeVarList1, treat.varname, outcome.varname, original=0, numKnot) {
  tryCatch ( 
    {
      
      data=NULL
      if (original==1){
        data=dataOR
      } else {
        
        treatID=which(dataOR[, treat.varname]==1)
        controlID=which(dataOR[, treat.varname]==0)
        data = dataOR[c( sample(treatID,replace=T), sample(controlID,replace=T) ),]  ##in the simulation studies, we did stratified bootstraps, results were similar with random bootstraps 
        #data = dataOR[sample(1:dim(dataInput)[1],replace=T),]  ###random bootstraps  ##
      }
      
      
      treatInd = dataOR[, treat.varname]  ###indicator of treated subjects in the original dataset
      treatBoot = data[, treat.varname]  ###indicator of treated subjects in the bootstrap sample
      Yobs = dataOR[, outcome.varname]
      
      ######## propensity score model ########################################################### 
      propen.model=NULL
      propen.model=formulaF(varList=propenVarList, y.name=treat.varname)
      
      
      ######## propensity score model ########################################################### 
      model2a=glm(propen.model, data=data, family="binomial", control = list(maxit = 50))
      prob.boot=predict(model2a, newdata=data, type="response")
      temp2a=NULL
      temp2a=as.numeric(treatBoot==1)*prob.boot + as.numeric(treatBoot==0)*(1-prob.boot)  ###estimate propensity of observed treatment
      data$pslogit = log(temp2a/(1-temp2a))   ####spline on the logit of propensity
      
      
      ####fit the prediction model for all the models
      ###for the outcome model separate model y0 and y1
      ##############################################################
      ##############################################################
      ###use equally spaced fixed knots assuming K knots
      data0=data[data$A0==0,]
      pspp0=pencompFit(y.varname=outcome.varname, x.varnames=outcomeVarList0, propen.score=data0$pslogit, data=data0, num.knot=numKnot)
      
      ###imputing the missing potential outcomes
      newData0=dataOR
      newData0[, treat.varname]=0
      predict2a=1 - predict(model2a, newdata=newData0, type="response")
      newData0$pslogit=log(predict2a/(1-predict2a))
      
      imputed0 = imputeF(newdata=newData0, model=pspp0, x.varnames=outcomeVarList0, propen.score.new=newData0$pslogit)
      mean(imputed0)
      
      ##########################################################################################################################
      ###prediction model for Y under treatment using treated data
      data1=data[data$A0==1,]
      pspp1=pencompFit(y.varname=outcome.varname, x.varnames=outcomeVarList1, propen.score=data1$pslogit, data=data1, num.knot=numKnot)
      
      ###imputing the missing potential outcomes
      newData1=dataOR
      newData1[, treat.varname]=1
      predict2a=predict(model2a, newdata=newData1, type="response")
      newData1$pslogit=log(predict2a/(1-predict2a))
      
      imputed1 = imputeF(newdata=newData1, model=pspp1, x.varnames=outcomeVarList1, propen.score.new=newData1$pslogit)
      mean(imputed1)
      
      ##############keep the observed outcome the same#####################
      imputed1[treatInd==1]=Yobs[treatInd==1]
      imputed0[treatInd==0]=Yobs[treatInd==0]
      
      #######include everyone
      ATE=c(mean(imputed1-imputed0), var(imputed1-imputed0)/nrow(dataOR) )
      
      
      return( ATE ) ###estimate and variance of ATE
    }, error=function(e) return(NA) )
}



