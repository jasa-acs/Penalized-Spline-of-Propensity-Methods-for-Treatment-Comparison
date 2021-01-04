##################################################################################################################
##################################################################################################################
### linear outcome simulation scenario
### after finishing running the simulations, we combine the results to obtain estimates for bias,
#relative bias, sd, RMSE, 95% coverage, avgerage confidence Width, number of simulation, mean bootstrap standard error.

rm(list=ls())

source("codes/twoTimePointSimulation/FiguresandTables/processResults.R")

sampleSizeNum=c(200, 500, 1000)

coef1Var=c(-0.5, -0.8, -0.8)
coef2Var=c(-0.1, -0.1, -0.5)
coef3Var=c(0.2, 0.6, 1.1)

####this script produces a file containing correct, misPred, misWeight for each simulation scenario
for (methodFolder in c("/PENCOMP_Results/", "/AIPTW_Results/", "/gcompute_Results/", "/IPTW_Results/")) {
  
for(outcome in c("LinearOutcome", "NonLinearOutcome")){
  
  
  if(outcome=="LinearOutcome"){
    
    truth11=22.35
    truth10=11.17
    truth01=10.45
    
  } else {
    
    truth11=25.31
    truth10=12.69
    truth01=10.57
    
  }
  
  
  
  for (numSize in 1:length(sampleSizeNum)) {
    
    sampleSize=sampleSizeNum[numSize]
    
    DIREC=paste0("M:/Private/Rpackage_v3/resubmission_round4/codes/twoTimePointSimulation/sampleSize", sampleSize, "/", outcome, methodFolder) ###directory where simulation results are stored
    
    for(var2 in 1:length(coef1Var))  {
      
      coef1=coef1Var[var2]
      coef2=coef2Var[var2]
      coef3=coef3Var[var2]
      
      
      ############################################################################################################
      ################ both models are correct #######################################################################################
      fileNames=list.files(DIREC, pattern=paste("correct_coef1_", coef1, "coef2_", coef2, "coef3_", coef3, ".txt", sep=""))
      
      psppMixed=NULL
      for(k in 1:length(fileNames)) {
        tempResult=read.table(paste(DIREC, fileNames[k], sep=""),header=T, sep="\t")
        tempResult=tempResult[which(tempResult[,1] != 0),]
        psppMixed=rbind(psppMixed,tempResult)
      }
      
      correctFile=psppMixed
      correct = processTotal(resultTotal=correctFile, truth11=truth11, truth10=truth10, truth01=truth01)
      
      
      ######################################### prediction models are incorrect #################################################################
      ###############################################################################################################################
      ############################################################################################################
      if(methodFolder=="/IPTW_Results/"){
        misPred=NULL
      } else {
      fileNames=list.files(DIREC,  pattern=paste("misPred_coef1_", coef1, "coef2_", coef2, "coef3_", coef3, ".txt", sep=""))
      
      psppMixed=NULL
      for(k in 1:length(fileNames)) {
        tempResult=read.table(paste(DIREC, fileNames[k], sep=""),header=T, sep="\t")
        tempResult=tempResult[which(tempResult[,1] != 0),]
        psppMixed=rbind(psppMixed,tempResult)
      }
      
      result=NULL
      result=psppMixed
      
      misPredFile=psppMixed
      misPred = processTotal(resultTotal=misPredFile, truth11=truth11, truth10=truth10, truth01=truth01)
      }
      
      ######################################### propensity models are incorrect #################################################################
      ###############################################################################################################################
      #####################################################################################################################################
      if(methodFolder=="/gcompute_Results/"){
        misWeight=NULL
      } else {
      fileNames=list.files(DIREC, pattern=paste("misWeight_coef1_", coef1, "coef2_", coef2, "coef3_", coef3, ".txt", sep=""))
      
      psppMixed=NULL
      for(k in 1:length(fileNames)) {
        tempResult=read.table(paste(DIREC, fileNames[k], sep=""),header=T, sep="\t")
        tempResult=tempResult[which(tempResult[,1] != 0),]
        psppMixed=rbind(psppMixed,tempResult)
      }
      
      misWeightFile=psppMixed
      misWeight = processTotal(resultTotal=misWeightFile, truth11=truth11, truth10=truth10, truth01=truth01)
      }
      
      tempOut=rbind(correct, misPred, misWeight) 
      
      colnames(tempOut)=rep(c("bias", "rel bias", "sd", "RMSE", "95% coverage","avg Width", "numSim", "meanBoot"), 3)
      options(digits=2)
      
      write.table(tempOut, paste(DIREC,  "resultOutput_", "coef1_", coef1, "coef2_", coef2, "coef3_", coef3, ".txt",sep=""), sep="\t", quote=F, row.names=T, col.names=T)
      
      }
    
    }
  }
}

