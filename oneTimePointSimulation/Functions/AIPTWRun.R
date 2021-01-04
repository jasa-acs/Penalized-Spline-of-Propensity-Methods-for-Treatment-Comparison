#################################################################################################
#################################################################################################
#### repeat the simulations for a large number of times
####change the parameters for different simulation scenarios as described in the paper 

####this will calculate AIPTW estimate, standard error and 95% confidence interval for x number of simulated datasets

###numSim: number of simulation
###numT: number of complete datasets used for estimating standard error and confidence interval
###linear: true or false ( for linear or nonlinear in covariates, specifying the two cases I looked at in the simulation studies)
###modelType: correct, misPred, misWeight, for the three cases considered in the study
###treat.varname: name of treatment indicator
###outcome.varname: name of outcome variable
##numKnot: number of knots in forming the linear basis functions
###numSim-number of simulations, numT-number of bootstrap samples to obtain standard error and confidence interval


rm(list=ls())
library(nlme)
library("mgcv")
library("MASS")
library("rootSolve") 

###step 1: load all the functions

DIREC="codes/oneTimePointSimulation/"

###functions for simulating data and additional functions used for processing outputs
###load the functions
source(paste(DIREC, "Functions/simulateData.R", sep=""))
source(paste(DIREC, "Functions/modelSpec.R", sep="")) 
source(paste(DIREC, "Functions/truth.R", sep="")) 
source(paste(DIREC, "Functions/formulaConstruct.R", sep=""))  
source(paste(DIREC, "Functions/addFun.R", sep=""))   

###functions for AIPTW
source(paste(DIREC, "Functions/AIPTW.R", sep=""))   ###AIPTW estimate for one dataset
source(paste(DIREC, "Functions/gcomputeFunc.R", sep="")) ###additional g function for computing g computation estimates



treat.varname="A0"
outcome.varname="y"


numSim=100   ###takes integer values-number of simulated datasets
numT=500     ##takes integer values  -number of bootstrap samples for each dataset
sampleSize=200   #, 500, 1000 
linear="true"   ##"false" for linear vs nonlinear outcome models 
modelType="correct" #; "misPred", "misWeight" for model specification

level="high"  ###"moderate", or "low" for degree of confounding
baselineVar=c("L1", "L2", "L3")  ###baseline covariates
numCarlo=2000  ###number of monte carlo simulations used in estimating the augmentation term


outcomeModel=NULL
if(linear=="true"){
  outcomeModel="/LinearOutcome/"
} else {
  outcomeModel="/NonLinearOutcome/"
}

DIRECOUT=paste0(DIREC, "sampleSize", sampleSize, outcomeModel)

tableAIPTW=matrix(NA, nrow=numSim, ncol=4)  ###estimate, sd and 95% confidence interval lower and upper bounds



######################## note for level
if(level=="high"){
  
  beta1=1.5
  beta2=1.5
  beta3=0.75
  
} else if(level=="moderate"){
  
  beta1=1
  beta2=1
  beta3=0.5    
  
} else {
  
  beta1=0.1
  beta2=0.1
  beta3=0.05
  
}


######################################

for(k in 1:numSim)
{
  tryCatch (
    {
      
      print(k)   
      
      ###simulate the dataset
      simdat=simulateData(sampleSize=sampleSize, level, seed.num=k, linear=linear)
      
      modOut=modelSpec(modelType=modelType, linear=linear)
      outcomeVarList0=modOut$modely0
      outcomeVarList1=modOut$modely1
      propenVarList=modOut$propen.model
      
      ###########################without reestimation ################
      ###
      aiptwEst=NULL
      aiptwEst=AIPTW(dataInput=simdat, baselineVar=baselineVar, linear=linear, propenVarList=propenVarList, outcomeVarList0=outcomeVarList0, outcomeVarList1=outcomeVarList1,
                     treat.varname=treat.varname, outcome.varname=outcome.varname,
                     original=1, numCarlo=numCarlo)
      
      ########bootstrap samples to estimate the se and confidence intervals##########
      mulResult=NULL
      mulResult=replicate(numT, AIPTW(dataInput=simdat, baselineVar=baselineVar, linear=linear, propenVarList=propenVarList, outcomeVarList0=outcomeVarList0, outcomeVarList1=outcomeVarList1,
                                      treat.varname=treat.varname, outcome.varname=outcome.varname,
                                      original=0, numCarlo=numCarlo) )
      
      
      tableAIPTW[k, ]=c( aiptwEst[1], sd(mulResult, na.rm = T), aiptwEst[1] + c(-1,1) * 1.96 * sd(mulResult, na.rm = T) )  
      
      write.table(tableAIPTW, paste(DIRECOUT, "AIPTW_Results/", modelType, "_", "beta1_", beta1, "beta2_", beta2, "beta3_", beta3, ".txt",sep=""),
                  row.name=F, quote=F,
                  col.names = c("estimate", "se", "lower","upper"),
                  sep="\t")
      
      
    }
    ,
    error=function(e) { }
  )
}




