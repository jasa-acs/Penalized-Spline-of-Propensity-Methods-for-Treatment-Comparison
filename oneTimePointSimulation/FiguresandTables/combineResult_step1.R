##################################################################################################################
##################################################################################################################
######## 

### after finishing running the simulations, we combine the results to obtain estimates for bias,
#relative bias, sd, RMSE, 95% coverage, avgerage confidence Width, number of simulation, mean bootstrap standard error.
###for each of the methods-pencomp, iptw, gcomputation and aiptw, run the following scripts

##########here I store the output files in the same folders as the simulation results 

rm(list=ls())
source("codes/oneTimePointSimulation/FiguresandTables/processResults.R")

sampleSizeNum=c(200, 500, 1000)


beta1Var=c(1.5, 1, 0.1)
beta2Var=c(1.5, 1, 0.1)
beta3Var=c(0.75, 0.5, 0.05)


#################################combine pencomp and AIPTW results ########################

for (methodFolder in c("/PENCOMP_Results/", "/AIPTW_Results/", "/gcompute_Results/", "/IPTW_Results/")) {

for(outcome in c("LinearOutcome", "NonLinearOutcome")){
  
  if(outcome=="LinearOutcome"){
    truth11=5
  } else {
    truth11=9
  } 

  
####this script produces a file containing correct, misPred, misWeight for each simulation scenario
  for (numSize in 1:length(sampleSizeNum)) {

    sampleSize=sampleSizeNum[numSize]
    DIREC=paste0("M:/Private/Rpackage_v3/resubmission_round4/codes/oneTimePointSimulation/sampleSize", sampleSize, "/", outcome, methodFolder) ###directory where simulation results are stored

    for(var2 in 1:length(beta1Var))  {
      
      beta1=beta1Var[var2]
      beta2=beta2Var[var2]
      beta3=beta3Var[var2]
      
      ############################################################################################################
      ################ both propensity and prediction models are correct #######################################################################################
      ###combine all the simulation together
      fileNames=list.files(DIREC, pattern=paste("correct_beta1_", beta1, "beta2_", beta2, "beta3_", beta3, ".txt", sep=""))
      
      ##combine all the simulations together
      psppMixed=NULL
      for(k in 1:length(fileNames)) {
        tempResult=read.table(paste(DIREC, fileNames[k], sep=""),header=T, sep="\t")
        tempResult=tempResult[which(tempResult[,1] != 0),]
        psppMixed=rbind(psppMixed,tempResult)
      }
      
      correctFile=psppMixed
      correct = processResult(result=correctFile, truth=truth11)

      
      ######################################### prediction model is incorrect #################################################################
      ###############################################################################################################################
      if(methodFolder=="/IPTW_Results/"){
        misPred=NULL
      } else {
      
      fileNames=list.files(DIREC, pattern=paste("misPred_beta1_", beta1, "beta2_", beta2, "beta3_", beta3, ".txt", sep=""))
      
      psppMixed=NULL
      for(k in 1:length(fileNames)) {
        tempResult=read.table(paste(DIREC, fileNames[k], sep=""),header=T, sep="\t")
        tempResult=tempResult[which(tempResult[,1] != 0),]
        psppMixed=rbind(psppMixed,tempResult)
      }

      misPredFile=psppMixed
      misPred = processResult(result=misPredFile, truth=truth11)
      
      } 

      ######################################### propensity score model is incorrect #################################################################
      ###############################################################################################################################
      ############################################################################################################
      if(methodFolder=="/gcompute_Results/"){
        misWeight=NULL
      } else {
        
      fileNames=list.files(DIREC, pattern=paste("misWeight_beta1_", beta1, "beta2_", beta2, "beta3_", beta3, ".txt", sep=""))
      
      psppMixed=NULL
      for(k in 1:length(fileNames)) {
        tempResult=read.table(paste(DIREC, fileNames[k], sep=""),header=T, sep="\t")
        tempResult=tempResult[which(tempResult[,1] != 0),]
        psppMixed=rbind(psppMixed,tempResult)
      }
      
      misWeightFile=psppMixed
      misWeight = processResult(result=misWeightFile, truth=truth11)
      
      } 
      
      tempOut=rbind(correct, misPred, misWeight) 
      
      colnames(tempOut)=c(rep(c("bias", "rel bias", "sd", "RMSE", "95% coverage","avg Width"), 1),"numSim", "meanBoot11")
      
      options(digits=2)
      
      write.table(tempOut, paste(DIREC, "resultOutput_", "beta1_", beta1, "beta2_", beta2, "beta3_", beta3, ".txt",sep=""), sep="\t", quote=F, row.names=T, col.names=T)
      
    
      }
  
    }
  }
}






