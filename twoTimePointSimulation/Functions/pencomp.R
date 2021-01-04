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
###treat.varnameA0--first treatment variable name
####treat.varnameA1--second treatment variable name
###outcome.varname--final outcome of interest
###original=1 for original dataset, 0 for bootstrap sample
###numCarlo=2000: number of runs in g computation 
##numKnot-number of knots in forming basis functions

pencomp=function(dataInput, linear, propenVarListA0.Init, propenVarListA1.Init,
                 outcomeVarListL1.Init, outcomeVarListL1b.Init,
                 outcomeVarListY11.Init, outcomeVarListY10.Init,
                 outcomeVarListY01.Init, outcomeVarListY00.Init,
                 intermediate.varname1, intermediate.varname2,
                 treat.varnameA0, treat.varnameA1, outcome.varname, original=1, numKnot) {
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
      

      ######## propensity models ########################################################### 
      model2aF=formulaF(varList=propenVarListA0.Init, y.name=treat.varnameA0)
      model2a=glm(model2aF, data=data, family="binomial") 
      temp2a=treatIndA0 * model2a$fitted.values + (1-treatIndA0) * (1-model2a$fitted.values)
      
      data$L1L0 = data$L1-data$L0
      data$L1bL0b = data$L1b-data$L0b
      model2bF = formulaF(varList=propenVarListA1.Init, y.name=treat.varnameA1)
      model2b = glm(model2bF, data=data, family="binomial") ### correctly specified treatment mechanism 
      temp2b = treatIndA1 * model2b$fitted.values + (1-treatIndA1) * (1-model2b$fitted.values)
      
      demo=temp2a * temp2b 
      data$sw = log(demo/(1-demo)) ###in logit scale, 2nd time point
      data$pslogit = log(temp2a/(1-temp2a))   ####in logit scale, 1st time point
      
      ##############################################################
      ###use equally spaced fixed knots assuming K knots
      data0=data[treatIndA0==0,]
      pspp0=pencompFit(y.varname=intermediate.varname1, x.varnames=outcomeVarListL1.Init[1], propen.score=data0$pslogit, data=data0, num.knot=numKnot)
      pspp0b=pencompFit(y.varname=intermediate.varname2, x.varnames=outcomeVarListL1b.Init, propen.score=data0$pslogit, data=data0, num.knot=numKnot)
      
      
      ############# IMPUTE THE MISSING VALUES IN THE ORIGINAL DATA SET#####################################
      newData0=data.frame(L0=data$L0, L0b=data$L0b, L1=rep(NA, dim(data)[1]), L1b=rep(NA, dim(data)[1]), A0=rep(0, dim(data)[1]))
      newData0$id=1
      predict2a=1-predict(model2a, newData0, type="response")
      newData0$pslogit=log(predict2a/(1-predict2a))  ###in logit scale
      
      newData0$L1=imputeF(newdata=newData0, model=pspp0, x.varnames=outcomeVarListL1.Init[1], propen.score.new=newData0$pslogit)
      newData0$L1[which(treatIndA0==0)]=data$L1[which(treatIndA0==0)]
      newData0$L1b=imputeF(newdata=newData0, model=pspp0b, x.varnames=outcomeVarListL1b.Init, propen.score.new=newData0$pslogit)
      newData0$L1b[which(treatIndA0==0)]=data$L1b[which(treatIndA0==0)]
      
      
      
      ##############################################################
      ###use equally spaced fixed knots assuming K knots
      data1=data[treatIndA0==1,]
      pspp1=pencompFit(y.varname=intermediate.varname1, x.varnames=outcomeVarListL1.Init[1], propen.score=data1$pslogit, data=data1, num.knot=numKnot)
      pspp1b=pencompFit(y.varname=intermediate.varname2, x.varnames=outcomeVarListL1b.Init, propen.score=data1$pslogit, data=data1, num.knot=numKnot)
      
      
      ############# IMPUTE THE MISSING VALUES IN THE ORIGINAL DATA SET#####################################
      newData1=data.frame(L0=data$L0, L0b=data$L0b, L1=rep(NA, dim(data)[1]), L1b=rep(NA, dim(data)[1]), A0=rep(0, dim(data)[1]))
      newData1$id=1
      predict2a=predict(model2a, newData1, type="response")
      newData1$pslogit=log(predict2a/(1-predict2a))  ###in logit scale
      
      newData1$L1=imputeF(newdata=newData1, model=pspp1, x.varnames=outcomeVarListL1.Init[1], propen.score.new=newData1$pslogit)
      newData1$L1[which(treatIndA0==1)]=data$L1[which(treatIndA0==1)]
      newData1$L1b=imputeF(newdata=newData1, model=pspp1b, x.varnames=outcomeVarListL1b.Init, propen.score.new=newData1$pslogit)
      newData1$L1b[which(treatIndA0==1)]=data$L1b[which(treatIndA0==1)]
      
      
      ####################IMPUTING THE MISSING VALUES IN THE ORIGINAL DATA SET################################
      newData00=data.frame(L0=data$L0, L0b=data$L0b, L1=newData0$L1, L1b=newData0$L1b, A0=rep(0, dim(data)[1]), A1=rep(0, dim(data)[1]))
      newData00$L1L0=newData00$L1-newData00$L0
      newData00$L1bL0b=newData00$L1b-newData00$L0b
      predict2a=1-predict(model2a, newData00, type="response")
      predict2b=1-predict(model2b, newData00, type="response")
      newData00$sw=log((predict2a * predict2b)/(1-predict2a * predict2b))  ###in logit scale
      newData00$id=1
      newData00$L1L1bInt = newData00$L1 * newData00$L1b ###interaction terms between L1 and L1b
      
      
      newData10=data.frame(L0=data$L0, L0b=data$L0b, L1=newData1$L1, L1b=newData1$L1b, A0=rep(1, dim(data)[1]), A1=rep(0, dim(data)[1]))
      newData10$L1L0=newData10$L1-newData10$L0
      newData10$L1bL0b=newData10$L1b-newData10$L0b
      predict2a=predict(model2a, newData10, type="response")    
      predict2b=1-predict(model2b, newData10, type="response")
      newData10$sw=log((predict2a * predict2b)/(1-predict2a * predict2b))  ###in logit scale
      newData10$id=1
      newData10$L1L1bInt = newData10$L1 * newData10$L1b ###interaction terms between L1 and L1b
      
      newData01=data.frame(L0=data$L0, L0b=data$L0b, L1=newData0$L1, L1b=newData0$L1b, A0=rep(0, dim(data)[1]), A1=rep(1, dim(data)[1]))
      newData01$L1L0=newData01$L1-newData01$L0
      newData01$L1bL0b=newData01$L1b-newData01$L0b
      predict2a=1-predict(model2a, newData01, type="response")    
      predict2b=predict(model2b, newData01, type="response")
      newData01$sw=log((predict2a * predict2b)/(1-predict2a * predict2b))  ###in logit scale
      newData01$id=1
      newData01$L1L1bInt = newData01$L1 * newData01$L1b ###interaction terms between L1 and L1b
      
      newData11=data.frame(L0=data$L0, L0b=data$L0b, L1=newData1$L1, L1b=newData1$L1b, A0=rep(1, dim(data)[1]), A1=rep(1, dim(data)[1]))
      newData11$L1L0=newData11$L1-newData11$L0
      newData11$L1bL0b=newData11$L1b-newData11$L0b
      predict2a=predict(model2a, newData11, type="response")    
      predict2b=predict(model2b, newData11, type="response")
      newData11$sw=log((predict2a * predict2b)/(1-predict2a * predict2b))  ###in logit scale
      newData11$id=1
      newData11$L1L1bInt = newData11$L1 * newData11$L1b ###interaction terms between L1 and L1b
      
      
      #######################################################################################################
      ###use equally spaced fixed knots assuming K knots
      data00=data[treatIndA0==0 & treatIndA1==0,]
      data00$L1L1bInt = data00$L1 * data00$L1b ###interaction terms between L1 and L1b
      pspp00=pencompFit(y.varname=outcome.varname, x.varnames=outcomeVarListY00.Init, propen.score=data00$sw, data=data00, num.knot=numKnot)
      
      
      ############################################################################################
      ###use equally spaced fixed knots assuming K knots
      data01=data[treatIndA0==0 & treatIndA1==1,]
      data01$L1L1bInt = data01$L1 * data01$L1b ###interaction terms between L1 and L1b
      pspp01=pencompFit(y.varname=outcome.varname, x.varnames=outcomeVarListY01.Init, propen.score=data01$sw, data=data01, num.knot=numKnot)
      
      
      #########################################################################################
      ###use equally spaced fixed knots assuming K knots
      data10=data[treatIndA0==1 & treatIndA1==0,]
      data10$L1L1bInt = data10$L1 * data10$L1b ###interaction terms between L1 and L1b
      pspp10=pencompFit(y.varname=outcome.varname, x.varnames=outcomeVarListY01.Init, propen.score=data10$sw, data=data10, num.knot=numKnot)
      
      
      #########################################################################################
      data11=data[treatIndA0==1 & treatIndA1==1,]
      data11$L1L1bInt = data11$L1 * data11$L1b ###interaction terms between L1 and L1b
      pspp11=pencompFit(y.varname=outcome.varname, x.varnames=outcomeVarListY01.Init, propen.score=data11$sw, data=data11, num.knot=numKnot)
      
      
      ########imputing the missing potential outcomes
      potential00=imputeF(newdata=newData00, model=pspp00, x.varnames=outcomeVarListY00.Init, propen.score.new=newData00$sw)
      potential00[which(treatIndA0==0 & treatIndA1==0)]=data[which(treatIndA0==0 & treatIndA1==0), outcome.varname]
      
      potential01=imputeF(newdata=newData01, model=pspp01, x.varnames=outcomeVarListY01.Init, propen.score.new=newData01$sw)
      potential01[which(treatIndA0==0 & treatIndA1==1)]=data[which(treatIndA0==0 & treatIndA1==1), outcome.varname]
      
      potential10=imputeF(newdata=newData10, model=pspp10, x.varnames=outcomeVarListY10.Init, propen.score.new=newData10$sw)
      potential10[which(treatIndA0==1 & treatIndA1==0)]=data[which(treatIndA0==1 & treatIndA1==0), outcome.varname]
      
      potential11=imputeF(newdata=newData11, model=pspp11, x.varnames=outcomeVarListY11.Init, propen.score.new=newData11$sw)
      potential11[which(treatIndA0==1 & treatIndA1==1)]=data[which(treatIndA0==1 & treatIndA1==1), outcome.varname]
      
      estimate11=potential11-potential00
      estimate10=potential10-potential00
      estimate01=potential01-potential00
      
      return( c(mean(estimate11), var(estimate11)/length(estimate11), 
                mean(estimate10), var(estimate10)/length(estimate10),
                mean(estimate01), var(estimate01)/length(estimate01)) )
      
    }, error=function(e) return( NA ) )

}


