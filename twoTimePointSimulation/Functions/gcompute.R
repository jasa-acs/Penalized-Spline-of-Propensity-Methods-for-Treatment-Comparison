#####this calculate g computation estimate
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
###numCarlo=2000: number of runs in g computation 

gcompute= function(dataInput, propenVarListA0.Init, propenVarListA1.Init,
                   
                   outcomeVarListL1.Init, outcomeVarListL1b.Init,
                   outcomeVarListY11.Init, outcomeVarListY10.Init,
                   outcomeVarListY01.Init, outcomeVarListY00.Init,
                   
                   intermediate.varname1, intermediate.varname2,
                   treat.varnameA0, treat.varnameA1, outcome.varname, original=1, numCarlo=2000, baselineVar) { 
  
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
      
      
      ###########first intermediate outcome
      ###specify the models for intermediate outcomes
      model_0 = model_1 = formulaF(varList=outcomeVarListL1.Init, y.name=intermediate.varname1)
      mod2_0=lm(model_0, data=data[treatIndA0==0,]) 
      mod2_1=lm(model_1, data=data[treatIndA0==1,]) 
      
      model_0 = model_1 = formulaF(varList=outcomeVarListL1b.Init, y.name=intermediate.varname2)
      mod2_0b=lm(model_0, data=data[treatIndA0==0,]) 
      mod2_1b=lm(model_1, data=data[treatIndA0==1,]) 
      
      ####specify the outcome model
      model_00 = formulaF(varList=outcomeVarListY00.Init, y.name=outcome.varname)
      model_01 = formulaF(varList=outcomeVarListY01.Init, y.name=outcome.varname)
      model_10 = formulaF(varList=outcomeVarListY10.Init, y.name=outcome.varname)
      model_11 = formulaF(varList=outcomeVarListY11.Init, y.name=outcome.varname)
      
      mod3_00=lm(model_00, data=data[treatIndA0==0 & treatIndA1==0, ]) 
      mod3_01=lm(model_01, data=data[treatIndA0==0 & treatIndA1==1, ]) 
      mod3_10=lm(model_10, data=data[treatIndA0==1 & treatIndA1==0, ]) 
      mod3_11=lm(model_11, data=data[treatIndA0==1 & treatIndA1==1, ]) 
      
      
      ############################################################### 
      ###simulate potential outcomes under each treatment regime
      Gcomp11=gcomputeFunc(inter.mod1=mod2_1, inter.mod2=mod2_1b, final.mod=mod3_11, data=data, baselineVar=baselineVar, numCarlo=numCarlo,
                           treat.varnameA0=treat.varnameA0, treat.varnameA1=treat.varnameA1,
                           treatSeq=c(1,1) )
      
      
      
      ############################################################### 
      ###simulate potential outcomes under each treatment regime
      Gcomp10=gcomputeFunc(inter.mod1=mod2_1, inter.mod2=mod2_1b, final.mod=mod3_10, data=data, baselineVar=baselineVar, numCarlo=numCarlo,
                           treat.varnameA0=treat.varnameA0, treat.varnameA1=treat.varnameA1,
                           treatSeq=c(1,0) )
      
      
      ############################################################### 
      ###simulate potential outcomes under each treatment regime
      Gcomp01=gcomputeFunc(inter.mod1=mod2_0, inter.mod2=mod2_0b, final.mod=mod3_01, data=data, baselineVar=baselineVar, numCarlo=numCarlo,
                           treat.varnameA0=treat.varnameA0, treat.varnameA1=treat.varnameA1,
                           treatSeq=c(0,1) )
      
      ############################################################### 
      ###simulate potential outcomes under each treatment regime
      Gcomp00=gcomputeFunc(inter.mod1=mod2_0, inter.mod2=mod2_0b, final.mod=mod3_00, data=data, baselineVar=baselineVar, numCarlo=numCarlo,
                           treat.varnameA0=treat.varnameA0, treat.varnameA1=treat.varnameA1,
                           treatSeq=c(0,0) )
      
      estFinal_msm11=mean(Gcomp11[, outcome.varname])-mean(Gcomp00[, outcome.varname])  ###estimate 
      estFinal_msm10=mean(Gcomp10[, outcome.varname])-mean(Gcomp00[, outcome.varname])  #
      estFinal_msm01=mean(Gcomp01[, outcome.varname])-mean(Gcomp00[, outcome.varname])  #
      
      
      return( c(estFinal_msm11, estFinal_msm10, estFinal_msm01) )
      
      
    }, error=function(e) return(NA) )
  
} 

