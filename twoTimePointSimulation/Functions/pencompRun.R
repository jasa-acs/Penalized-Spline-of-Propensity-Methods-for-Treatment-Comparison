#################################################################################################
#################################################################################################
#### repeat the simulations for a large number of times
####change the parameters for different simulation scenarios as described in the paper 
####this will calculate pencomp estimate, standard error and 95% confidence interval for x number of simulated datasets


rm(list=ls())
library(nlme)
library("mgcv")
library("MASS")

###step 1: load all the functions

DIREC="codes/twoTimePointSimulation/"  ###root directory where functions and results would be stored
###functions for simulating data and additional functions used for processing outputs
funcLoc="Functions/"
source(paste(DIREC, funcLoc, "simulateData.R", sep=""))  ###simualte dataset
source(paste(DIREC, funcLoc, "modelSpec.R", sep=""))  ##model specification
source(paste(DIREC, funcLoc, "truth.R", sep=""))    ###estimating the true treatment effects 
source(paste(DIREC, funcLoc, "formulaConstruct.R", sep="")) ###construct formua
source(paste(DIREC, funcLoc, "addFun.R", sep=""))   ###additional functions

source(paste(DIREC, funcLoc, "pencompFit.R", sep=""))  ###fit penalized spline model with truncated linear basis 
source(paste(DIREC, funcLoc, "pencomp.R", sep=""))  ###pencomp estimate for one dataset


treat.varnameA0="A0"
treat.varnameA1="A1"
outcome.varname="y"
intermediate.varname1="L1"
intermediate.varname2="L1b"


numKnot=15  ###number of knots in the penalized spline
numSim=1000 ##number of simulated datasets
numT=200 ##number of complete datasets created
sampleSize=200 ##, 500, 1000
linear="false"  ###linear vs nonlinear
modelType="correct" ##misPred, and misWeight
level="high"

outcomeModel=NULL
if(linear=="true"){
  outcomeModel="/LinearOutcome/"
} else{
  outcomeModel="/NonLinearOutcome/"
}

DIRECOUT=paste0(DIREC, "sampleSize", sampleSize, outcomeModel)

tablePENCOMP=matrix(NA, nrow=numSim, ncol=4*3)  ###Rubin's combining rules four entries, estimate, se and confidence intervals (lower and upper) and 3 estimands
###delta11=E(y11-y00), delta10=E(y10-y00) and delta01=E(y01-y00)

############################################
coef1=coef2=coef3=NULL

if(level=="high"){
  
  coef1=-0.8
  coef2=-0.5
  coef3=1.1
  
} else if(level=="moderate"){
  
  coef1=-0.8
  coef2=-0.1
  coef3=0.6
  
} else {
  
  coef1=-0.5
  coef2=-0.1
  coef3=0.2
  
}

##############################################

for(d in 1:numSim)
{
  tryCatch ( 
    {
      
      print(d)   
      
      set.seed(d)
      ###simulate the dataset
      simdat=simulateData(sampleSize=sampleSize, level=level, seed.num=d, linear=linear)
      
      
      ##propensity model
      modOut=modelSpec(modelType=modelType, linear=linear)
      propenVarListA0.Init=modOut$propenVarListA0.Init
      propenVarListA1.Init=modOut$propenVarListA1.Init
      
      ###models for intermediate outcomes
      outcomeVarListL1.Init=modOut$outcomeVarListL1.Init
      outcomeVarListL1b.Init=modOut$outcomeVarListL1b.Init
      
      ###models for final outcome
      outcomeVarListY11.Init=modOut$outcomeVarListY11.Init
      outcomeVarListY10.Init=modOut$outcomeVarListY10.Init
      outcomeVarListY01.Init=modOut$outcomeVarListY01.Init
      outcomeVarListY00.Init=modOut$outcomeVarListY00.Init
      
      #######bootstrap samples to get standard errors
      ########bootstrap samples
      mulResult=NULL
      mulResult=replicate(numT, pencomp(dataInput=simdat, linear=linear, propenVarListA0.Init=propenVarListA0.Init,
                                        propenVarListA1.Init=propenVarListA1.Init,
                                        outcomeVarListL1.Init=outcomeVarListL1.Init, 
                                        outcomeVarListL1b.Init=outcomeVarListL1b.Init,
                                        outcomeVarListY11.Init=outcomeVarListY11.Init,
                                        outcomeVarListY10.Init=outcomeVarListY10.Init,
                                        outcomeVarListY01.Init=outcomeVarListY01.Init,
                                        outcomeVarListY00.Init=outcomeVarListY00.Init,
                                        intermediate.varname1=intermediate.varname1,
                                        intermediate.varname2=intermediate.varname2,
                                        treat.varnameA0=treat.varnameA0,
                                        treat.varnameA1=treat.varnameA1, outcome.varname=outcome.varname, original=0, numKnot=numKnot) )
      
      
      
      for(ind in c(0, 1, 2) ){
        tablePENCOMP[d, (1:4)+(4*ind)] =processPENCOMP(mulResult[(1:2)+(2*ind),]) 
      }
      
      
      
      write.table(tablePENCOMP, paste(DIRECOUT, "PENCOMP_Results/", modelType, "_", "coef1_", coef1, "coef2_", coef2, "coef3_", coef3,  ".txt",sep=""), row.name=F, quote=F,
                  col.names = rep(c("estimate", "se", "lower","upper"), 3),
                  sep="\t")
      
      
    }, error=function(e) { } )              
}




