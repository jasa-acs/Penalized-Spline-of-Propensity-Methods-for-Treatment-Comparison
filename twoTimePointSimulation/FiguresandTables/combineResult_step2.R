
##################################################################################################################
##################################################################################################################
###### COMBINE RESULTS FROM DIFFERENT CONFOUNDING##########


###after combining the simulation results for each simulation scenario using PENCOMPResult.step1.R
## This script aggregates results so there is only one file for each outcome model and samplesize (for example labeled as sampleSize200_LinearOutcome.txt)

###for each of the methods-pencomp, iptw, g computation and aiptw, run the following script

rm(list=ls())

sampleSizeNum=c(200, 500, 1000)

DIREC="codes/twoTimePointSimulation/" ###where the files are stored
DIRECOUT="codes/twoTimePointSimulation/FiguresandTables"  ####where the output files are stored

for(method in c("/PENCOMP_Results/", "/IPTW_Results/", "/gcompute_Results/", "/AIPTW_Results/")) {
  
for (outcome in c("NonLinearOutcome", "LinearOutcome") ){
    
  for (numSize in 1:length(sampleSizeNum)) {
  
  
    sampleSize=sampleSizeNum[numSize]
    
    
    low=read.table(paste(DIREC, "sampleSize", sampleSize,"/", outcome,  method,"resultOutput_coef1_-0.5coef2_-0.1coef3_0.2.txt", sep=""), header=T, sep="\t")
    moderate=read.table(paste(DIREC, "sampleSize", sampleSize,"/", outcome,  method,"resultOutput_coef1_-0.8coef2_-0.1coef3_0.6.txt", sep=""), header=T, sep="\t")
    high=read.table(paste(DIREC, "sampleSize", sampleSize,"/", outcome,  method,"resultOutput_coef1_-0.8coef2_-0.5coef3_1.1.txt", sep=""), header=T, sep="\t")
    
    tempOut=NULL
    
    if(method == "/gcompute_Results/"){
      
      tempOut=rbind(low[which(row.names(low)=="correct"),], moderate[which(row.names(moderate)=="correct"),], high[which(row.names(high)=="correct"),],
                    low[which(row.names(low)=="misPred"),], moderate[which(row.names(moderate)=="misPred"),], high[which(row.names(high)=="misPred"),]) 
      
      colnames(tempOut)=rep(c("bias", "rel bias", "sd", "RMSE", "95% coverage","avg Width", "numSim", "meanBoot"), 3)
      
      row.names(tempOut)=c("lowCorrect", "modCorrect", "highCorrect", "lowMisPred", "modMispred", "highMispred")
      
    } else if(method == "/IPTW_Results/"){
      
      tempOut=rbind(low[which(row.names(low)=="correct"),], moderate[which(row.names(moderate)=="correct"),], high[which(row.names(high)=="correct"),],
                    low[which(row.names(low)=="misWeight"),], moderate[which(row.names(moderate)=="misWeight"),], high[which(row.names(high)=="misWeight"),]) 
      
      colnames(tempOut)=rep(c("bias", "rel bias", "sd", "RMSE", "95% coverage","avg Width", "numSim", "meanBoot"), 3)
      
      row.names(tempOut)=c("lowCorrect", "modCorrect", "highCorrect", "lowmisWeight", "modmisWeight", "highmisWeight")
      
      
    } else {
      
      tempOut=rbind(low[which(row.names(low)=="correct"),], moderate[which(row.names(moderate)=="correct"),], high[which(row.names(high)=="correct"),],
                    low[which(row.names(low)=="misPred"),], moderate[which(row.names(moderate)=="misPred"),], high[which(row.names(high)=="misPred"),],
                    low[which(row.names(low)=="misWeight"),], moderate[which(row.names(moderate)=="misWeight"),], high[which(row.names(high)=="misWeight"),]) 
      
      
      colnames(tempOut)=rep(c("bias", "rel bias", "sd", "RMSE", "95% coverage","avg Width", "numSim", "meanBoot"), 3)
      
      row.names(tempOut)=c("lowCorrect", "modCorrect", "highCorrect", "lowMisPred", "modMispred", "highMispred",
                           "lowMisWeight", "modMisWeight", "highMisWeight")
 
    }
  
    options(digits=2)
    
    write.table(tempOut, paste(DIRECOUT, method, "sampleSize", sampleSize, "_", outcome, ".txt", sep=""), sep="\t", quote=F, row.names=T, col.names=T)
    }
  }
}

