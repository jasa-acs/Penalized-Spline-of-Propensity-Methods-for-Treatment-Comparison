#################################################################################################
#################################################################################################
#### repeat the simulations for a large number of times
####change the parameters for different simulation scenarios as described in the paper 

###numSim: number of simulation
###numT: number of complete datasets used for estimating standard error and confidence interval
###linear: true or false ( for linear or nonlinear in covariates, specifying the two cases I looked at in the simulation studies)
###modelType: correct, misPred, misWeight, for the three cases considered in the study
###treat.varname: name of treatment indicator
###outcome.varname: name of outcome variable

rm(list=ls())
library(nlme)
library("mgcv")
library("MASS")

###step 1: load all the functions

DIREC="codes/oneTimePointSimulation/"  ###directory where files and functions are stored
###functions for simulating data and additional functions used for processing outputs
source(paste(DIREC, "Functions/simulateData.R", sep=""))
source(paste(DIREC, "Functions/modelSpec.R", sep="")) 
source(paste(DIREC, "Functions/truth.R", sep="")) 
source(paste(DIREC, "Functions/formulaConstruct.R", sep=""))  
source(paste(DIREC, "Functions/addFun.R", sep=""))   

####functions for IPTW
source(paste(DIREC, "Functions/IPTW.R", sep=""))  ## IPTW estimate for one dataset

treat.varname="A0"
outcome.varname="y"


numSim=200   ###takes integer values-number of simulated datasets
numT=200     ##takes integer values  -number of bootstrap samples for each dataset
sampleSize=200   #, 500, 1000 
linear="true"   ##"false" for linear vs nonlinear outcome models 
modelType="correct" #; "misPred", "misWeight" for model specification

level="high"  ###"moderate", or "low" for degree of confounding
baselineVar=c("L1", "L2", "L3")  ###baseline covariates


outcomeModel=NULL
if(linear=="true"){
  outcomeModel="/LinearOutcome/"
} else {
  outcomeModel="/NonLinearOutcome/"
}

DIRECOUT=paste0(DIREC, "sampleSize", sampleSize, outcomeModel)

tableIPTW=matrix(NA, nrow=numSim, ncol=4)  ###four entries, estimate, se and confidence intervals (lower and upper)



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

############################################

for(d in 1:numSim)
{
  tryCatch (
    {
      
      print(d)   
      
      set.seed(d)
      ###simulate the dataset
      simdat=simulateData(sampleSize=sampleSize, level=level, seed.num=d, linear=linear)
      
      modOut=modelSpec(modelType=modelType, linear=linear)
      propenVarList=modOut$propen.model
      
      ###########################without reestimation ################
      ###
      iptwEst=NULL
      iptwEst=IPTW(dataInput=simdat, propenVarList=propenVarList, treat.varname=treat.varname, outcome.varname=outcome.varname, original=1)
      
      ########bootstrap samples to estimate se and confidence intervals
      mulResult=NULL
      mulResult=replicate(numT, IPTW(dataInput=simdat, propenVarList=propenVarList, treat.varname=treat.varname, outcome.varname=outcome.varname, original=0) )
      
      
      tableIPTW[d, ]=c( iptwEst, sd(mulResult, na.rm = T), iptwEst + c(-1,1)*1.96*sd(mulResult, na.rm = T) )  
      
      
      write.table(tableIPTW, paste(DIRECOUT, "IPTW_Results/", modelType, "_", "beta1_", beta1, "beta2_", beta2, "beta3_", beta3, ".txt",sep=""), row.name=F, quote=F,
                  col.names = c("estimate", "se", "lower","upper"),
                  sep="\t")
      
      
    }
    ,
    error=function(e) { }
  )
}




