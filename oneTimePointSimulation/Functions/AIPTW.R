### function to calculate the doubly robust coefficients from the estimating equations 
###dataInput-the dataset
###propenVarlist-the variables in the propensity model
###treat.varname-names of the treatment 
###outcome.varname-name of the outcome variable
##outcomeVarList0-variables included in the outcome model for Y0
##outcomeVarList1-varialbes included in the outcome model for Y1
##original-1 for calculating treatment effect estimate on the original dataset, instead of bootstrap samples (for calculating standard error)
###numCarlo-number of monte carlo runs in g computation and estimating the augmentation term in AIPTW
###baseList-variables at baseline, draw from baseline empirical distributions (in the simulations there are 3 baseline covariates for one time point treatment)


AIPTW = function(dataInput, baselineVar, linear, propenVarList, outcomeVarList0, outcomeVarList1, treat.varname, outcome.varname, original, numCarlo) { 
  
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
      
      
      ######## MSM weight########################################################### 
      ###propensity model on numerator
      propen.model1a=formulaF(varList=c(1), y.name=treat.varname)
      model1a=glm(propen.model1a, data=data, family="binomial") 
      temp1a=treatInd * model1a$fitted.values + (1-treatInd) * (1-model1a$fitted.values)
      
      ###propensity model on denominator
      propen.model2a=formulaF(varList=propenVarList, y.name=treat.varname)
      model2a=glm(propen.model2a, data=data, family="binomial") 
      temp2a=treatInd * model2a$fitted.values + (1-treatInd) * (1-model2a$fitted.values)
      
      weight=temp1a/temp2a 
      
      summary(weight)
      
      #### use g-computation to compute beta(Qn) 
      ###prediction model for Y1 and Y0
      modely0 = formulaF(varList=outcomeVarList0, y.name=outcome.varname)
      modely1 = formulaF(varList=outcomeVarList1, y.name=outcome.varname)
      
      mod3_1=lm(modely1, data=data[treatInd==1,]) 
      mod3_0=lm(modely0, data=data[treatInd==0,]) 
      
      ## numCarlo number of runs in G-computation 
      ## numCarlo number of runs in G-computation 
      ############################################################### 
      ###simulate potential outcomes under each treatment A0=1
      Gcomp1=gcomputeFunc(final.mod=mod3_1, data=data, baselineVar=baselineVar, numCarlo=numCarlo, treat.varname=treat.varname, outcome.varname=outcome.varname, treat=1)
      
      
      ############################################################### 
      ###simulate potential outcomes under each treatment A0=0
      Gcomp0=gcomputeFunc(final.mod=mod3_0, data=data, baselineVar=baselineVar, numCarlo=numCarlo, treat.varname=treat.varname, outcome.varname=outcome.varname, treat=0)
      
      Gcomp=rbind(Gcomp1, Gcomp0) 
      dim(Gcomp) 
      
      ##get the treatment coefficients in the MSM model #### 
      designMat=cbind(Gcomp[, which(names(Gcomp)==treat.varname)])
      designY=Gcomp[, which(names(Gcomp)==outcome.varname)]
      msmModel=lm(designY ~ designMat, data=Gcomp) 
      summary(msmModel)
      
      
      d=matrix(0, nrow=dim(dataInput)[1], ncol=2) 
      for(u in 1:dim(dataInput)[1]) { 
        
        currentData=NULL ##store baseline coariate
        ### generate L0 from empirical distribution of L0 obtained from the data
        for(ind in 1:length(baselineVar)){
          draw_g=NULL
          draw_g=cbind(draw_g, rep(as.numeric(data[u, baselineVar[ind] ]), numCarlo))
          colnames(draw_g)=baselineVar[ind]
          currentData=cbind(currentData, draw_g)
        }
        
        currentData=data.frame(currentData, A0_g=data[u, treat.varname])
        names(currentData)[which(names(currentData) == "A0_g")]=treat.varname
        
        ####squared terms ###############
        currentData[ , "L1_sq"] = currentData[, "L1"]^2
        currentData[ , "L2_sq"] = currentData[, "L2"]^2
        currentData[ , "L3_sq"] = currentData[, "L3"]^2
        currentData[, "L1L2"] = currentData[, "L1"] * currentData[, "L2"]
        
        ####potential outcome under control
        y_mean0 = as.vector(predict(mod3_0, newdata=currentData))
        y_draw0=rnorm(numCarlo, mean=y_mean0, sd=summary(mod3_0)$sigma)
        
        ###potential outcome under treatment
        y_mean1 = as.vector(predict(mod3_1, newdata=currentData))
        y_draw1=rnorm(numCarlo, mean=y_mean1, sd=summary(mod3_1)$sigma)
        
        ####
        firstTreat = currentData[, treat.varname]
        y_draw=y_draw0 * (1 - firstTreat) + y_draw1 * firstTreat
        
        currentData[, outcome.varname]=y_draw
        
        temp1a_mc_=predict(model1a, currentData, type="response") 
        temp1a_mc=firstTreat * temp1a_mc_ + (1 - firstTreat) * (1-temp1a_mc_)
        
        temp2a_mc_=predict(model2a, currentData, type="response") 
        temp2a_mc=firstTreat * temp2a_mc_ + (1 - firstTreat)* (1-temp2a_mc_)
        
        weight_mc=temp1a_mc/temp2a_mc 
        
        temp=weight_mc * ( y_draw - (cbind(rep(1, numCarlo), firstTreat) %*% msmModel$coef) ) 
        d[u,]=c( mean(temp), mean(firstTreat * temp) ) 
        
      } 
      
      ########################################################################################################### 
      ### first expectation given history up to time poin 0 (baseline) 
      d_=matrix(0, nrow=dim(dataInput)[1], ncol=2) 
      
      for (u in 1:dim(dataInput)[1]) { 
        
        currentData=NULL ##store baseline coariate
        ### generate L0 from empirical distribution of L0 obtained from the data
        for(ind in 1:length(baselineVar)){
          draw_g=NULL
          draw_g=cbind(draw_g, rep(as.numeric(data[u, baselineVar[ind] ]), numCarlo))
          colnames(draw_g)=baselineVar[ind]
          currentData=cbind(currentData, draw_g)
        }
        
        currentData=data.frame(currentData, A0_g=data[u, treat.varname])
        names(currentData)[which(names(currentData) == "A0_g")]=treat.varname
        
        ####squared terms ###############
        currentData[ , "L1_sq"] = currentData[, "L1"]^2
        currentData[ , "L2_sq"] = currentData[, "L2"]^2
        currentData[ , "L3_sq"] = currentData[, "L3"]^2
        currentData[, "L1L2"] = currentData[, "L1"] * currentData[, "L2"]
        
        prob = predict(model2a, currentData, type = "response")  ###probability of treatment at first time point
        
        A0_draw=rbinom(numCarlo, 1, prob ) ### draw treatment at time point 2 
        currentData[, treat.varname]=A0_draw
        
        ####potential outcome under control
        y_mean0 = as.vector(predict(mod3_0, newdata=currentData))
        y_draw0=rnorm(numCarlo, mean=y_mean0, sd=summary(mod3_0)$sigma)
        
        ###potential outcome under treatment
        y_mean1 = as.vector(predict(mod3_1, newdata=currentData))
        y_draw1=rnorm(numCarlo, mean=y_mean1, sd=summary(mod3_1)$sigma)
        
        ####
        firstTreat = currentData[, treat.varname]
        y_draw=y_draw0 * (1 - firstTreat) + y_draw1 * firstTreat
        
        currentData[, outcome.varname]=y_draw
        
        temp1a_mc_=predict(model1a, currentData, type="response") 
        temp1a_mc=firstTreat * temp1a_mc_ + (1 - firstTreat) * (1-temp1a_mc_)
        
        temp2a_mc_=predict(model2a, currentData, type="response") 
        temp2a_mc=firstTreat * temp2a_mc_ + (1 - firstTreat)* (1-temp2a_mc_)
        
        weight_mc=temp1a_mc/temp2a_mc 
        
        temp=weight_mc * ( y_draw - (cbind(rep(1, numCarlo), firstTreat) %*% msmModel$coef) ) 
        d_[u,]=c( mean(temp), mean(firstTreat * temp) ) 
        
      } 
      
      ######### Newton-Raphson to solve the estimating equations 
      #install.packages("rootSolve") 
      #library("rootSolve") 
      
      firstTreat = NULL
      firstTreat = data[, treat.varname]
      m=cbind(rep(1,sampleSize), firstTreat)
      m2= t(d-d_)%*% rep(1,sampleSize)
      
      D<-function(beta){
        m<-cbind( rep(1,sampleSize), firstTreat ) %*% beta
        return(
          t(cbind(1, firstTreat)) %*% (weight*(data[, outcome.varname] - m))-m2
        )
      }
      
      oldmybeta<-c(0,5)
      
      repeat{
        mybeta <- oldmybeta - solve( gradient( D, oldmybeta) ) %*% D(oldmybeta)
        if (max(abs(oldmybeta-mybeta))<1e-8) break
        oldmybeta<-mybeta
      }
      
      result=c(1,1)%*%mybeta-c(1,0)%*%mybeta
      
      return(result) 
      
    }, error=function(e) return(NA) )
  
} 

