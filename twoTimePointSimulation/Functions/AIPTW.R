###
#####this calculate aiptw estimate
####dataInput--input dataset
###linear--linear or nonlinear in covariates
##propenVarListA0.Init--propensity model at first time point
####propenVarListA1.Init--propensity model at second time point
####outcomeVarListL1.Init--intermediate outcome model for L1
###outcomeVarListL1b.Init--intermediate outcome model for L1b
###outcomeVarListY11.Init--final outcome model for Y under treatment sequence of 11
###outcomeVarListY10.Init--final outcome model for Y under treatment sequence of 10
####outcomeVarListY01.Init--final outcome model for Y under treatment sequence of 01
####outcomeVarListY00.Init--final outcome model for Y under treatment sequence of 00
###intermediate.varname1--name of first intermediate outcome L1 in the simulation study
###intermediate.varname2--name of second intermediate outcome L1b
###baselineVar--names of baseline covariates
###treat.varnameA0--first treatment variable name
####treat.varnameA1--second treatment variable name
###outcome.varname--final outcome of interest
###original=1 for original dataset, 0 for bootstrap sample
###numCarlo=2000: number of replicates to calculate the augmentation term

AIPTW= function(dataInput, linear, propenVarListA0.Init, propenVarListA1.Init,
                outcomeVarListL1.Init, outcomeVarListL1b.Init,
                outcomeVarListY11.Init, outcomeVarListY10.Init,
                outcomeVarListY01.Init, outcomeVarListY00.Init,
                intermediate.varname1, intermediate.varname2, baselineVar, 
                treat.varnameA0, treat.varnameA1, outcome.varname, original=1, numCarlo=2000) { 
  
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
      #########calculate stabilized weights##########################################
      ############### numerator
      model1aF=formulaF(varList=c(1), y.name=treat.varnameA0)
      model1a=glm(model1aF, data=data, family="binomial") 
      temp1a=treatIndA0 * model1a$fitted.values + (1-treatIndA0) * (1-model1a$fitted.values)
      
      model1bF=formulaF(varList=c(treat.varnameA0), y.name=treat.varnameA1)
      model1b=glm(model1bF, data=data, family="binomial") 
      temp1b=treatIndA1 * model1b$fitted.values + (1-treatIndA1) * (1-model1b$fitted.values)
      
      numer=temp1a * temp1b 
      
      ################ denominator 
      model2aF=formulaF(varList=propenVarListA0.Init, y.name=treat.varnameA0)
      model2a=glm(model2aF, data=data, family="binomial") 
      temp2a=treatIndA0 * model2a$fitted.values + (1-treatIndA0) * (1-model2a$fitted.values)
      
      data$L1L0 = data$L1-data$L0
      data$L1bL0b = data$L1b-data$L0b
      model2bF = formulaF(varList=propenVarListA1.Init, y.name=treat.varnameA1)
      model2b = glm(model2bF, data=data, family="binomial") 
      temp2b = treatIndA1 * model2b$fitted.values + (1-treatIndA1) * (1-model2b$fitted.values)
      
      demo=temp2a * temp2b 
      weight=numer/demo  
    
      
      ###########first intermediate outcome
      ###specify the models for intermediate outcomes
      model_0 = model_1 = formulaF(varList=outcomeVarListL1.Init, y.name=intermediate.varname1)
      mod2_0=lm(model_0, data=data[treatIndA0==0,]) ###model for intermediate outcome L1 under first treatment A0=0
      mod2_1=lm(model_1, data=data[treatIndA0==1,]) ### under A0=1
      
      model_0 = model_1 = formulaF(varList=outcomeVarListL1b.Init, y.name=intermediate.varname2)
      mod2_0b=lm(model_0, data=data[treatIndA0==0,]) ###model for intermediate outcome L1b under first treatment A0=0
      mod2_1b=lm(model_1, data=data[treatIndA0==1,])  ###under A0=1
      
      ####specify the outcome model
      model_00 = formulaF(varList=outcomeVarListY00.Init, y.name=outcome.varname)
      model_01 = formulaF(varList=outcomeVarListY01.Init, y.name=outcome.varname)
      model_10 = formulaF(varList=outcomeVarListY10.Init, y.name=outcome.varname)
      model_11 = formulaF(varList=outcomeVarListY11.Init, y.name=outcome.varname)
      
      mod3_00=lm(model_00, data=data[treatIndA0==0 & treatIndA1==0, ])   ###final outcome model under treatment sequence of 11
      mod3_01=lm(model_01, data=data[treatIndA0==0 & treatIndA1==1, ])   ###under treatment sequence of 01
      mod3_10=lm(model_10, data=data[treatIndA0==1 & treatIndA1==0, ])   ##under treatment sequence of 10
      mod3_11=lm(model_11, data=data[treatIndA0==1 & treatIndA1==1, ])  ###under treatment sequence of 11
      
      
      ############################################################### 
      ###simulate potential outcomes under each treatment regime
      Gcomp11=gcomputeFunc(inter.mod1=mod2_1, inter.mod2=mod2_1b, final.mod=mod3_11, data=data, baselineVar=baselineVar, numCarlo=10000,
                           treat.varnameA0=treat.varnameA0, treat.varnameA1=treat.varnameA1,
                           treatSeq=c(1,1) )
      
      
      
      ############################################################### 
      ###simulate potential outcomes under each treatment regime
      Gcomp10=gcomputeFunc(inter.mod1=mod2_1, inter.mod2=mod2_1b, final.mod=mod3_10, data=data, baselineVar=baselineVar, numCarlo=10000,
                           treat.varnameA0=treat.varnameA0, treat.varnameA1=treat.varnameA1,
                           treatSeq=c(1,0) )
      
      
      ############################################################### 
      ###simulate potential outcomes under each treatment regime
      Gcomp01=gcomputeFunc(inter.mod1=mod2_0, inter.mod2=mod2_0b, final.mod=mod3_01, data=data, baselineVar=baselineVar, numCarlo=10000,
                           treat.varnameA0=treat.varnameA0, treat.varnameA1=treat.varnameA1,
                           treatSeq=c(0,1) )
      
      ############################################################### 
      ###simulate potential outcomes under each treatment regime
      Gcomp00=gcomputeFunc(inter.mod1=mod2_0, inter.mod2=mod2_0b, final.mod=mod3_00, data=data, baselineVar=baselineVar, numCarlo=10000,
                           treat.varnameA0=treat.varnameA0, treat.varnameA1=treat.varnameA1,
                           treatSeq=c(0,0) )
      
      Gcomp=rbind(Gcomp11, Gcomp10, Gcomp01, Gcomp00) 
      dim(Gcomp) 
      
      ##get the treatment coefficients in the MSM model #### 
      designMat=cbind(Gcomp[, which(names(Gcomp)==treat.varnameA0)], Gcomp[, which(names(Gcomp)==treat.varnameA1)],
                      Gcomp[, which(names(Gcomp)==treat.varnameA0)] * Gcomp[, which(names(Gcomp)==treat.varnameA1)])
      designY=Gcomp[, which(names(Gcomp)==outcome.varname)]
      msmModel=lm(designY ~ designMat, data=Gcomp) 
      summary(msmModel)
      
      
      
      ####################################################################  
      ########## first term in the augmentation term #####################
      ########here I calculate the augmentation term, in a two time point treatment, there are 4 sub parts in the augmentation term
      d=matrix(NA, nrow = nrow(data), ncol=4)
      for (numSubj in 1:nrow(data)){  ###for each individual, I simulate numCarlo replicates
        
        currentData=NULL ##store baseline coariate
        ### generate L0 from empirical distribution of L0 obtained from the data
        for(ind in 1:length(baselineVar)){
          draw_g=NULL
          draw_g=cbind(draw_g, rep(as.numeric(data[numSubj, baselineVar[ind] ]), numCarlo))
          colnames(draw_g)=baselineVar[ind]
          currentData=cbind(currentData, draw_g)
          
        }
        currentData=data.frame(currentData, A0_g=data[numSubj, treat.varnameA0])
        names(currentData)[which(names(currentData) == "A0_g")]=treat.varnameA0
        
        
        ### intermediate outcome L1 first expectation given history up to time point 0 (baseline) 
        L1_mean0=predict(mod2_0, newdata = currentData) 
        L1_draw0=rnorm(numCarlo, mean=L1_mean0, sd=summary(mod2_0)$sigma)   
        
        L1_mean1=predict(mod2_1, newdata = currentData) 
        L1_draw1=rnorm(numCarlo, mean=L1_mean1, sd=summary(mod2_1)$sigma)  
        
        firstTreat = currentData[, treat.varnameA0]
        L1_draw=L1_draw0*(1 - firstTreat) + L1_draw1 * firstTreat
        currentData[, intermediate.varname1] = L1_draw
        
        
        ### intermediate outcome L1b first expectation given history up to time point 0 (baseline) 
        L1b_mean0=predict(mod2_0b, newdata = currentData)
        L1b_draw0=rnorm(numCarlo, mean=L1b_mean0, sd=summary(mod2_0b)$sigma)   
        
        L1b_mean1=predict(mod2_1b, newdata = currentData)
        L1b_draw1=rnorm(numCarlo, mean=L1b_mean1, sd=summary(mod2_1b)$sigma)   
        
        L1b_draw=L1b_draw0*(1-firstTreat) + L1b_draw1*firstTreat
        currentData[, intermediate.varname2] = L1b_draw
        
        currentData$L1L0 = currentData$L1 - currentData$L0
        currentData$L1bL0b = currentData$L1b - currentData$L0b
        currentData$L1L1bInt = currentData$L1 * currentData$L1b
        prob=predict(model2b, newdata = currentData, type="response") 
        
        A1_draw=rbinom(numCarlo,1,prob) ### draw treatment at time point 2 
        currentData[, treat.varnameA1] = A1_draw
        
        y_mean00=predict(mod3_00, newdata = currentData) 
        y_draw00=rnorm(numCarlo, mean=y_mean00, sd=summary(mod3_00)$sigma) 
        
        ##########################################################
        y_mean01=predict(mod3_01, newdata = currentData) 
        y_draw01=rnorm(numCarlo, mean=y_mean01, sd=summary(mod3_01)$sigma) 
        
        #############################################################
        y_mean10=predict(mod3_10, newdata = currentData) 
        y_draw10=rnorm(numCarlo, mean=y_mean10, sd=summary(mod3_10)$sigma) 
        
        #############################################################
        y_mean11=predict(mod3_11, newdata = currentData) 
        y_draw11=rnorm(numCarlo, mean=y_mean11, sd=summary(mod3_11)$sigma) 
        
        secondTreat = currentData[, treat.varnameA1]
        currentData[, outcome.varname]=y_draw11 * as.numeric(firstTreat==1 & secondTreat==1) + y_draw10*as.numeric(firstTreat==1 & secondTreat==0) +
          y_draw01*as.numeric(firstTreat==0 & secondTreat==1) + y_draw00*as.numeric(firstTreat==0 & secondTreat==0)
        
        temp1a_mc_=predict(model1a, currentData, type="response") 
        temp1a_mc=firstTreat * temp1a_mc_ + (1-firstTreat)*(1-temp1a_mc_)
        
        temp1b_mc_=predict(model1b, currentData, type="response") 
        temp1b_mc=secondTreat * temp1b_mc_ + (1-secondTreat)*(1-temp1b_mc_)
        
        temp2a_mc_=predict(model2a, currentData, type="response") 
        temp2a_mc=firstTreat * temp2a_mc_ + (1-firstTreat)*(1-temp2a_mc_)
        
        temp2b_mc_=predict(model2b, currentData, type="response") 
        temp2b_mc=secondTreat * temp2b_mc_ + (1-secondTreat)*(1-temp2b_mc_)
        
        weight_mc=(temp1a_mc * temp1b_mc)/(temp2a_mc * temp2b_mc) 
        temp=weight_mc * (currentData[, outcome.varname] - (cbind(rep(1,numCarlo), firstTreat, secondTreat, firstTreat * secondTreat) %*% msmModel$coef)) 
        d[numSubj,]=c(mean(temp), mean(firstTreat * temp), mean(secondTreat * temp), mean(firstTreat * secondTreat * temp)) 
        
      }
      
      
      
      ##################################################################################################### 
      ### second expectation given history up to time point 1  
      ###################################################################################
      ##################################################################################
      
      d2=matrix(NA, nrow=nrow(data), ncol=4)
      
      for(numSubj in 1:nrow(data)){
        
        currentData=NULL
        ########### generate second intermediate outcome L1b 
        includeVar = NULL
        includeVar = c(baselineVar, intermediate.varname1, intermediate.varname2, treat.varnameA0, treat.varnameA1)
        for(ind in 1:length(includeVar)){
          draw_g=rep(as.numeric(data[numSubj, includeVar[ind]]), numCarlo)
          currentData=cbind(currentData, draw_g)
        }
        currentData=data.frame(currentData)
        names(currentData)=includeVar
        currentData$L1L0 = currentData$L1 - currentData$L0
        currentData$L1bL0b = currentData$L1b - currentData$L0b
        currentData$L1L1bInt = currentData$L1 * currentData$L1b
        
        ###########################
        y_mean00=predict(mod3_00, newdata = currentData) 
        y_draw00=rnorm(numCarlo, mean=y_mean00, sd=summary(mod3_00)$sigma) 
        
        ###########################
        y_mean01=predict(mod3_01, newdata = currentData) 
        y_draw01=rnorm(numCarlo, mean=y_mean01, sd=summary(mod3_01)$sigma) 
        
        ###########################
        y_mean10=predict(mod3_10, newdata = currentData) 
        y_draw10=rnorm(numCarlo, mean=y_mean10, sd=summary(mod3_10)$sigma) 
        
        
        ###########################
        y_mean11=predict(mod3_11, newdata = currentData) 
        y_draw11=rnorm(numCarlo, mean=y_mean11, sd=summary(mod3_11)$sigma) 
        
        firstTreat = currentData[, treat.varnameA0]
        secondTreat = currentData[, treat.varnameA1]
        currentData[, outcome.varname]=y_draw11 * as.numeric(firstTreat==1 & secondTreat==1) + y_draw10*as.numeric(firstTreat==1 & secondTreat==0) +
          y_draw01*as.numeric(firstTreat==0 & secondTreat==1) + y_draw00*as.numeric(firstTreat==0 & secondTreat==0)
        
        temp1a_mc_=predict(model1a, currentData, type="response") 
        temp1a_mc=firstTreat * temp1a_mc_ + (1-firstTreat)*(1-temp1a_mc_)
        
        temp1b_mc_=predict(model1b, currentData, type="response") 
        temp1b_mc=secondTreat * temp1b_mc_ + (1-secondTreat)*(1-temp1b_mc_)
        
        temp2a_mc_=predict(model2a, currentData, type="response") 
        temp2a_mc=firstTreat * temp2a_mc_ + (1-firstTreat)*(1-temp2a_mc_)
        
        temp2b_mc_=predict(model2b, currentData, type="response") 
        temp2b_mc=secondTreat * temp2b_mc_ + (1-secondTreat)*(1-temp2b_mc_)
        
        weight_mc=(temp1a_mc * temp1b_mc)/(temp2a_mc * temp2b_mc) 
        temp=weight_mc * (currentData[, outcome.varname] - (cbind(rep(1,numCarlo), firstTreat, secondTreat, firstTreat * secondTreat) %*% msmModel$coef)) 
        d2[numSubj,]=c(mean(temp), mean(firstTreat * temp), mean(secondTreat * temp), mean(firstTreat * secondTreat * temp)) 
        
      }
      
      d_sum=d+d2 
      
      
      #######################################################################
      #######################################################################
      ########################################################################################################### 
      ### first expectation given history up to time poin 0 (baseline) 
      d_=matrix(NA, nrow=nrow(data), ncol=4)
      for(numSubj in 1:nrow(data)){
        
        currentData=NULL
        ########### generate second intermediate outcome L1b 
        includeVar = NULL
        includeVar = baselineVar
        for(ind in 1:length(includeVar)){
          draw_g=rep(as.numeric(data[numSubj, includeVar[ind]]), numCarlo)
          currentData=cbind(currentData, draw_g)
        }
        currentData=data.frame(currentData)
        names(currentData)=includeVar
        
        #######################################################
        ######first generate the first treatment ###############
        prob1=predict(model2a, newdata = currentData, type = "response")
        A0_draw=rbinom(numCarlo, 1, prob = prob1) ### draw treatment at time point 2 
        currentData[, treat.varnameA0]=A0_draw
        
        
        ### intermediate outcome L1 first expectation given history up to time point 0 (baseline) 
        L1_mean0=predict(mod2_0, newdata = currentData) 
        L1_draw0=rnorm(numCarlo, mean=L1_mean0, sd=summary(mod2_0)$sigma)   
        
        L1_mean1=predict(mod2_1, newdata = currentData) 
        L1_draw1=rnorm(numCarlo, mean=L1_mean1, sd=summary(mod2_1)$sigma)  
        
        firstTreat = currentData[, treat.varnameA0]
        L1_draw=L1_draw0*(1 - firstTreat) + L1_draw1 * firstTreat
        currentData[, intermediate.varname1] = L1_draw
        
        
        ### intermediate outcome L1b first expectation given history up to time point 0 (baseline) 
        L1b_mean0=predict(mod2_0b, newdata = currentData)
        L1b_draw0=rnorm(numCarlo, mean=L1b_mean0, sd=summary(mod2_0b)$sigma)   
        
        L1b_mean1=predict(mod2_1b, newdata = currentData)
        L1b_draw1=rnorm(numCarlo, mean=L1b_mean1, sd=summary(mod2_1b)$sigma)   
        
        L1b_draw=L1b_draw0*(1-firstTreat) + L1b_draw1*firstTreat
        currentData[, intermediate.varname2] = L1b_draw
        
        currentData$L1L0 = currentData$L1 - currentData$L0
        currentData$L1bL0b = currentData$L1b - currentData$L0b
        currentData$L1L1bInt = currentData$L1 * currentData$L1b
        prob=predict(model2b, newdata = currentData, type="response") 
        
        A1_draw=rbinom(numCarlo,1,prob) ### draw treatment at time point 2 
        currentData[, treat.varnameA1] = A1_draw
        
        y_mean00=predict(mod3_00, newdata = currentData) 
        y_draw00=rnorm(numCarlo, mean=y_mean00, sd=summary(mod3_00)$sigma) 
        
        ##########################################################
        y_mean01=predict(mod3_01, newdata = currentData) 
        y_draw01=rnorm(numCarlo, mean=y_mean01, sd=summary(mod3_01)$sigma) 
        
        #############################################################
        y_mean10=predict(mod3_10, newdata = currentData) 
        y_draw10=rnorm(numCarlo, mean=y_mean10, sd=summary(mod3_10)$sigma) 
        
        #############################################################
        y_mean11=predict(mod3_11, newdata = currentData) 
        y_draw11=rnorm(numCarlo, mean=y_mean11, sd=summary(mod3_11)$sigma) 
        
        secondTreat = currentData[, treat.varnameA1]
        currentData[, outcome.varname]=y_draw11 * as.numeric(firstTreat==1 & secondTreat==1) + y_draw10*as.numeric(firstTreat==1 & secondTreat==0) +
          y_draw01*as.numeric(firstTreat==0 & secondTreat==1) + y_draw00*as.numeric(firstTreat==0 & secondTreat==0)
        
        temp1a_mc_=predict(model1a, currentData, type="response") 
        temp1a_mc=firstTreat * temp1a_mc_ + (1-firstTreat)*(1-temp1a_mc_)
        
        temp1b_mc_=predict(model1b, currentData, type="response") 
        temp1b_mc=secondTreat * temp1b_mc_ + (1-secondTreat)*(1-temp1b_mc_)
        
        temp2a_mc_=predict(model2a, currentData, type="response") 
        temp2a_mc=firstTreat * temp2a_mc_ + (1-firstTreat)*(1-temp2a_mc_)
        
        temp2b_mc_=predict(model2b, currentData, type="response") 
        temp2b_mc=secondTreat * temp2b_mc_ + (1-secondTreat)*(1-temp2b_mc_)
        
        weight_mc=(temp1a_mc * temp1b_mc)/(temp2a_mc * temp2b_mc) 
        temp=weight_mc * (currentData[, outcome.varname] - (cbind(rep(1,numCarlo), firstTreat, secondTreat, firstTreat * secondTreat) %*% msmModel$coef)) 
        d_[numSubj,]=c(mean(temp), mean(firstTreat * temp), mean(secondTreat * temp), mean(firstTreat * secondTreat * temp)) 
        
      }
      
      ########################################################## 
      ### second expectation given history up to time point 1 
      d2_=matrix(NA, nrow=nrow(data), ncol=4)
      
      for(numSubj in 1:nrow(data)){
        
        currentData=NULL
        ########### replicate numCarlo copies of prior history for each subject up to including second time point covariate history
        ##### before treatment assignment
        includeVar = NULL
        includeVar = c(baselineVar, intermediate.varname1, intermediate.varname2, treat.varnameA0)
        for(ind in 1:length(includeVar)){
          draw_g=rep(as.numeric(data[numSubj, includeVar[ind]]), numCarlo)
          currentData=cbind(currentData, draw_g)
        }
        currentData=data.frame(currentData)
        names(currentData)=includeVar
        
        #### simulate second treatment assignment based on subject prior treatment and covariate histories
        currentData$L1L0 = currentData$L1 - currentData$L0
        currentData$L1bL0b = currentData$L1b - currentData$L0b
        currentData$L1L1bInt = currentData$L1 * currentData$L1b
        prob = predict(model2b, currentData, type = "response")  ###probability of treatment at second time point 
        
        A1_draw=rbinom(numCarlo, 1, prob = prob) ### draw treatment at time point 2 
        currentData[, treat.varnameA1]=A1_draw
        
        y_mean00=predict(mod3_00, newdata = currentData) 
        y_draw00=rnorm(numCarlo, mean=y_mean00, sd=summary(mod3_00)$sigma) 
        
        ##########################################################
        y_mean01=predict(mod3_01, newdata = currentData) 
        y_draw01=rnorm(numCarlo, mean=y_mean01, sd=summary(mod3_01)$sigma) 
        
        #############################################################
        y_mean10=predict(mod3_10, newdata = currentData) 
        y_draw10=rnorm(numCarlo, mean=y_mean10, sd=summary(mod3_10)$sigma) 
        
        #############################################################
        y_mean11=predict(mod3_11, newdata = currentData) 
        y_draw11=rnorm(numCarlo, mean=y_mean11, sd=summary(mod3_11)$sigma) 
        
        firstTreat = currentData[, treat.varnameA0]
        secondTreat = currentData[, treat.varnameA1]
        currentData[, outcome.varname]=y_draw11 * as.numeric(firstTreat==1 & secondTreat==1) + y_draw10*as.numeric(firstTreat==1 & secondTreat==0) +
          y_draw01*as.numeric(firstTreat==0 & secondTreat==1) + y_draw00*as.numeric(firstTreat==0 & secondTreat==0)
        
        temp1a_mc_=predict(model1a, currentData, type="response") 
        temp1a_mc=firstTreat * temp1a_mc_ + (1-firstTreat)*(1-temp1a_mc_)
        
        temp1b_mc_=predict(model1b, currentData, type="response") 
        temp1b_mc=secondTreat * temp1b_mc_ + (1-secondTreat)*(1-temp1b_mc_)
        
        temp2a_mc_=predict(model2a, currentData, type="response") 
        temp2a_mc=firstTreat * temp2a_mc_ + (1-firstTreat)*(1-temp2a_mc_)
        
        temp2b_mc_=predict(model2b, currentData, type="response") 
        temp2b_mc=secondTreat * temp2b_mc_ + (1-secondTreat)*(1-temp2b_mc_)
        
        weight_mc=(temp1a_mc * temp1b_mc)/(temp2a_mc * temp2b_mc) 
        temp=weight_mc * (currentData[, outcome.varname] - (cbind(rep(1,numCarlo), firstTreat, secondTreat, firstTreat * secondTreat) %*% msmModel$coef)) 
        d2_[numSubj,]=c(mean(temp), mean(firstTreat * temp), mean(secondTreat * temp), mean(firstTreat * secondTreat * temp)) 
        
      }
      
      d_sum2=d_+d2_ 
      
      ######### Newton-Raphson to solve the estimating equations 
      #install.packages("rootSolve") 
      #library("rootSolve") 
      #################### solving the estimating equation to obtain doubly robust estimate
      
      firstTreat = as.vector(data[, treat.varnameA0])
      secondTreat = as.vector(data[, treat.varnameA1])
      
      m2= t(d_sum - d_sum2)%*% rep(1, sampleSize)  ###this is the augmentation term
      
      D<-function(beta){
        m<-cbind(rep(1,sampleSize), firstTreat, secondTreat, firstTreat * secondTreat)%*%beta
        return(
          t(cbind(1, firstTreat, secondTreat, firstTreat * secondTreat) ) %*% (weight * (data[, outcome.varname] - m )) - m2
        )
      }
      
      oldmybeta<-c(5,0,0,0)
      
      repeat{
        mybeta <- oldmybeta - solve( gradient(D, oldmybeta) ) %*% D(oldmybeta)
        if ( max( abs(oldmybeta-mybeta) )<1e-8 ) break
        oldmybeta <- mybeta
      }
      
      result<-c(c(1,1,1,1)%*%mybeta-c(1,0,0,0)%*%mybeta,
                c(1,1,0,0)%*%mybeta-c(1,0,0,0)%*%mybeta,
                c(1,0,1,0)%*%mybeta-c(1,0,0,0)%*%mybeta)
      
      return(result) 
      
    }, error=function(e) return(NA) )
  
} 