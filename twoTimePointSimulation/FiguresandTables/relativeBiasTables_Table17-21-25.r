########created 5/27/2017
##################################### bias as a proportion of RMSE of correct IPTW ##########################################################
rm(list=ls())


DIREC="codes/twoTimePointSimulation/FiguresandTables/" ###where files are stored after combineResult_step2.R

measure1="bias" ##for delta 11
measure2="bias.1"  ##for delta10
measure3="bias.2"  ##for delta01


measure1a="RMSE"  ##for delta11
measure2a="RMSE.1"  ##for delta10
measure3a="RMSE.2"  ##for delta01

sampleNum=c(200, 500, 1000)

for( num in 1:3) {
  
  sampleSize=sampleNum[num]
  
  pencompL=read.table(paste(DIREC, "PENCOMP_Results/sampleSize", sampleSize, "_LinearOutcome.txt", sep=""), header=T, sep="\t")
  pencompN=read.table(paste(DIREC, "PENCOMP_Results/sampleSize", sampleSize, "_NonLinearOutcome.txt", sep=""), header=T, sep="\t")
  
  gcomputeL=read.table(paste(DIREC, "gcompute_Results/sampleSize", sampleSize, "_LinearOutcome.txt", sep=""), header=T, sep="\t")
  gcomputeN=read.table(paste(DIREC, "gcompute_Results/sampleSize", sampleSize, "_NonLinearOutcome.txt", sep=""), header=T, sep="\t")
  
  
  iptwL=read.table(paste(DIREC, "IPTW_Results/sampleSize", sampleSize, "_LinearOutcome.txt", sep=""), header=T, sep="\t")
  iptwN=read.table(paste(DIREC, "IPTW_Results/sampleSize", sampleSize, "_NonLinearOutcome.txt", sep=""), header=T, sep="\t")
  
  
  aiptwL=read.table(paste(DIREC, "AIPTW_Results/sampleSize", sampleSize, "_LinearOutcome.txt", sep=""), header=T, sep="\t")
  aiptwN=read.table(paste(DIREC, "AIPTW_Results/sampleSize", sampleSize, "_NonLinearOutcome.txt", sep=""), header=T, sep="\t")
  
  
  bench1=c(iptwL[1,which(names(iptwL)==measure1a)],iptwL[2,which(names(iptwL)==measure1a)],iptwL[3,which(names(iptwL)==measure1a)], 
           iptwL[1,which(names(iptwL)==measure2a)],iptwL[2,which(names(iptwL)==measure2a)],iptwL[3,which(names(iptwL)==measure2a)], 
           iptwL[1,which(names(iptwL)==measure3a)],iptwL[2,which(names(iptwL)==measure3a)],iptwL[3,which(names(iptwL)==measure3a)])  ###benchmarks RMSE of correct IPTW
  
  bench2=c(iptwN[1,which(names(iptwN)==measure1a)],iptwN[2,which(names(iptwN)==measure1a)],iptwN[3,which(names(iptwN)==measure1a)], 
           iptwN[1,which(names(iptwN)==measure2a)],iptwN[2,which(names(iptwN)==measure2a)],iptwN[3,which(names(iptwN)==measure2a)], 
           iptwN[1,which(names(iptwN)==measure3a)],iptwN[2,which(names(iptwN)==measure3a)],iptwN[3,which(names(iptwN)==measure3a)])  ###benchmarks RMSE of correct IPTW
  
  ####################
  iptwBias1=rbind(c(iptwL[1,which(names(iptwL)==measure1)],iptwL[2,which(names(iptwL)==measure1)],iptwL[3,which(names(iptwL)==measure1)],  ### linear delta11
                    iptwL[1,which(names(iptwL)==measure2)],iptwL[2,which(names(iptwL)==measure2)],iptwL[3,which(names(iptwL)==measure2)], ##linear delta10
                    iptwL[1,which(names(iptwL)==measure3)],iptwL[2,which(names(iptwL)==measure3)],iptwL[3,which(names(iptwL)==measure3)])/bench1,  ##linear delta01
                  
                  c(iptwL[4,which(names(iptwL)==measure1)],iptwL[5,which(names(iptwL)==measure1)],iptwL[6,which(names(iptwL)==measure1)],  ### linear delta11
                    iptwL[4,which(names(iptwL)==measure2)],iptwL[5,which(names(iptwL)==measure2)],iptwL[6,which(names(iptwL)==measure2)], ##linear delta10
                    iptwL[4,which(names(iptwL)==measure3)],iptwL[5,which(names(iptwL)==measure3)],iptwL[6,which(names(iptwL)==measure3)])/bench1)  ##linear delta01
  
  
  
  iptwBias2=rbind(c(iptwN[1,which(names(iptwN)==measure1)],iptwN[2,which(names(iptwN)==measure1)],iptwN[3,which(names(iptwN)==measure1)],  ### linear delta11
                    iptwN[1,which(names(iptwN)==measure2)],iptwN[2,which(names(iptwN)==measure2)],iptwN[3,which(names(iptwN)==measure2)], ##linear delta10
                    iptwN[1,which(names(iptwN)==measure3)],iptwN[2,which(names(iptwN)==measure3)],iptwN[3,which(names(iptwN)==measure3)])/bench2,  ##linear delta01
                  
                  c(iptwN[4,which(names(iptwN)==measure1)],iptwN[5,which(names(iptwN)==measure1)],iptwN[6,which(names(iptwN)==measure1)],  
                    iptwN[4,which(names(iptwN)==measure2)],iptwN[5,which(names(iptwN)==measure2)],iptwN[6,which(names(iptwN)==measure2)], 
                    iptwN[4,which(names(iptwN)==measure3)],iptwN[5,which(names(iptwN)==measure3)],iptwN[6,which(names(iptwN)==measure3)])/bench2) 
  
  
  
  
  pencompBias1=rbind(c(pencompL[1,which(names(pencompL)==measure1)],pencompL[2,which(names(pencompL)==measure1)],pencompL[3,which(names(pencompL)==measure1)],  ### linear delta11
                       pencompL[1,which(names(pencompL)==measure2)],pencompL[2,which(names(pencompL)==measure2)],pencompL[3,which(names(pencompL)==measure2)], ##linear delta10
                       pencompL[1,which(names(pencompL)==measure3)],pencompL[2,which(names(pencompL)==measure3)],pencompL[3,which(names(pencompL)==measure3)])/bench1,##linear delta01
                     
                     c(pencompL[4,which(names(pencompL)==measure1)],pencompL[5,which(names(pencompL)==measure1)],pencompL[6,which(names(pencompL)==measure1)],  ### linear delta11
                       pencompL[4,which(names(pencompL)==measure2)],pencompL[5,which(names(pencompL)==measure2)],pencompL[6,which(names(pencompL)==measure2)], ##linear delta10
                       pencompL[4,which(names(pencompL)==measure3)],pencompL[5,which(names(pencompL)==measure3)],pencompL[6,which(names(pencompL)==measure3)])/bench1,##linear delta01
                     
                     c(pencompL[7,which(names(pencompL)==measure1)],pencompL[8,which(names(pencompL)==measure1)],pencompL[9,which(names(pencompL)==measure1)],  ### linear delta11
                       pencompL[7,which(names(pencompL)==measure2)],pencompL[8,which(names(pencompL)==measure2)],pencompL[9,which(names(pencompL)==measure2)], ##linear delta10
                       pencompL[7,which(names(pencompL)==measure3)],pencompL[8,which(names(pencompL)==measure3)],pencompL[9,which(names(pencompL)==measure3)])/bench1)##linear delta01
  
  
  pencompBias2=rbind(c(pencompN[1,which(names(pencompN)==measure1)],pencompN[2,which(names(pencompN)==measure1)],pencompN[3,which(names(pencompN)==measure1)],  ### linear delta11
                       pencompN[1,which(names(pencompN)==measure2)],pencompN[2,which(names(pencompN)==measure2)],pencompN[3,which(names(pencompN)==measure2)], ##linear delta10
                       pencompN[1,which(names(pencompN)==measure3)],pencompN[2,which(names(pencompN)==measure3)],pencompN[3,which(names(pencompN)==measure3)])/bench2,##linear delta01
                     
                     c(pencompN[4,which(names(pencompN)==measure1)],pencompN[5,which(names(pencompN)==measure1)],pencompN[6,which(names(pencompN)==measure1)],  ### linear delta11
                       pencompN[4,which(names(pencompN)==measure2)],pencompN[5,which(names(pencompN)==measure2)],pencompN[6,which(names(pencompN)==measure2)], ##linear delta10
                       pencompN[4,which(names(pencompN)==measure3)],pencompN[5,which(names(pencompN)==measure3)],pencompN[6,which(names(pencompN)==measure3)])/bench2,##linear delta01
                     
                     c(pencompN[7,which(names(pencompN)==measure1)],pencompN[8,which(names(pencompN)==measure1)],pencompN[9,which(names(pencompN)==measure1)],  ### linear delta11
                       pencompN[7,which(names(pencompN)==measure2)],pencompN[8,which(names(pencompN)==measure2)],pencompN[9,which(names(pencompN)==measure2)], ##linear delta10
                       pencompN[7,which(names(pencompN)==measure3)],pencompN[8,which(names(pencompN)==measure3)],pencompN[9,which(names(pencompN)==measure3)])/bench2)##linear delta01
  
  
  
  ####################
  aiptwBias1=rbind(c(aiptwL[1,which(names(aiptwL)==measure1)],aiptwL[2,which(names(aiptwL)==measure1)],aiptwL[3,which(names(aiptwL)==measure1)],  ### linear delta11
                     aiptwL[1,which(names(aiptwL)==measure2)],aiptwL[2,which(names(aiptwL)==measure2)],aiptwL[3,which(names(aiptwL)==measure2)], ##linear delta10
                     aiptwL[1,which(names(aiptwL)==measure3)],aiptwL[2,which(names(aiptwL)==measure3)],aiptwL[3,which(names(aiptwL)==measure3)])/bench1,##linear delta01
                   
                   c(aiptwL[4,which(names(aiptwL)==measure1)],aiptwL[5,which(names(aiptwL)==measure1)],aiptwL[6,which(names(aiptwL)==measure1)],  ### linear delta11
                     aiptwL[4,which(names(aiptwL)==measure2)],aiptwL[5,which(names(aiptwL)==measure2)],aiptwL[6,which(names(aiptwL)==measure2)], ##linear delta10
                     aiptwL[4,which(names(aiptwL)==measure3)],aiptwL[5,which(names(aiptwL)==measure3)],aiptwL[6,which(names(aiptwL)==measure3)])/bench1,##linear delta01
                   
                   c(aiptwL[7,which(names(aiptwL)==measure1)],aiptwL[8,which(names(aiptwL)==measure1)],aiptwL[9,which(names(aiptwL)==measure1)],  ### linear delta11
                     aiptwL[7,which(names(aiptwL)==measure2)],aiptwL[8,which(names(aiptwL)==measure2)],aiptwL[9,which(names(aiptwL)==measure2)], ##linear delta10
                     aiptwL[7,which(names(aiptwL)==measure3)],aiptwL[8,which(names(aiptwL)==measure3)],aiptwL[9,which(names(aiptwL)==measure3)])/bench1)##linear delta01
  
  
  aiptwBias2=rbind(c(aiptwN[1,which(names(aiptwN)==measure1)],aiptwN[2,which(names(aiptwN)==measure1)],aiptwN[3,which(names(aiptwN)==measure1)],  ### linear delta11
                     aiptwN[1,which(names(aiptwN)==measure2)],aiptwN[2,which(names(aiptwN)==measure2)],aiptwN[3,which(names(aiptwN)==measure2)], ##linear delta10
                     aiptwN[1,which(names(aiptwN)==measure3)],aiptwN[2,which(names(aiptwN)==measure3)],aiptwN[3,which(names(aiptwN)==measure3)])/bench2,##linear delta01
                   
                   c(aiptwN[4,which(names(aiptwN)==measure1)],aiptwN[5,which(names(aiptwN)==measure1)],aiptwN[6,which(names(aiptwN)==measure1)],  ### linear delta11
                     aiptwN[4,which(names(aiptwN)==measure2)],aiptwN[5,which(names(aiptwN)==measure2)],aiptwN[6,which(names(aiptwN)==measure2)], ##linear delta10
                     aiptwN[4,which(names(aiptwN)==measure3)],aiptwN[5,which(names(aiptwN)==measure3)],aiptwN[6,which(names(aiptwN)==measure3)])/bench2,##linear delta01
                   
                   c(aiptwN[7,which(names(aiptwN)==measure1)],aiptwN[8,which(names(aiptwN)==measure1)],aiptwN[9,which(names(aiptwN)==measure1)],  ### linear delta11
                     aiptwN[7,which(names(aiptwN)==measure2)],aiptwN[8,which(names(aiptwN)==measure2)],aiptwN[9,which(names(aiptwN)==measure2)], ##linear delta10
                     aiptwN[7,which(names(aiptwN)==measure3)],aiptwN[8,which(names(aiptwN)==measure3)],aiptwN[9,which(names(aiptwN)==measure3)])/bench2)##linear delta01
  
  
  
  ####################
  gcomputeBias1=rbind(c(gcomputeL[1,which(names(gcomputeL)==measure1)],gcomputeL[2,which(names(gcomputeL)==measure1)],gcomputeL[3,which(names(gcomputeL)==measure1)],  ### linear delta11
                        gcomputeL[1,which(names(gcomputeL)==measure2)],gcomputeL[2,which(names(gcomputeL)==measure2)],gcomputeL[3,which(names(gcomputeL)==measure2)], ##linear delta10
                        gcomputeL[1,which(names(gcomputeL)==measure3)],gcomputeL[2,which(names(gcomputeL)==measure3)],gcomputeL[3,which(names(gcomputeL)==measure3)])/bench1,  ##linear delta01
                      
                      c(gcomputeL[4,which(names(gcomputeL)==measure1)],gcomputeL[5,which(names(gcomputeL)==measure1)],gcomputeL[6,which(names(gcomputeL)==measure1)],  ### linear delta11
                        gcomputeL[4,which(names(gcomputeL)==measure2)],gcomputeL[5,which(names(gcomputeL)==measure2)],gcomputeL[6,which(names(gcomputeL)==measure2)], ##linear delta10
                        gcomputeL[4,which(names(gcomputeL)==measure3)],gcomputeL[5,which(names(gcomputeL)==measure3)],gcomputeL[6,which(names(gcomputeL)==measure3)])/bench1)  ##linear delta01
  
  
  
  gcomputeBias2=rbind(c(gcomputeN[1,which(names(gcomputeN)==measure1)],gcomputeN[2,which(names(gcomputeN)==measure1)],gcomputeN[3,which(names(gcomputeN)==measure1)],  ### linear delta11
                        gcomputeN[1,which(names(gcomputeN)==measure2)],gcomputeN[2,which(names(gcomputeN)==measure2)],gcomputeN[3,which(names(gcomputeN)==measure2)], ##linear delta10
                        gcomputeN[1,which(names(gcomputeN)==measure3)],gcomputeN[2,which(names(gcomputeN)==measure3)],gcomputeN[3,which(names(gcomputeN)==measure3)])/bench2,  ##linear delta01
                      
                      c(gcomputeN[4,which(names(gcomputeN)==measure1)],gcomputeN[5,which(names(gcomputeN)==measure1)],gcomputeN[6,which(names(gcomputeN)==measure1)],  ### linear delta11
                        gcomputeN[4,which(names(gcomputeN)==measure2)],gcomputeN[5,which(names(gcomputeN)==measure2)],gcomputeN[6,which(names(gcomputeN)==measure2)], ##linear delta10
                        gcomputeN[4,which(names(gcomputeN)==measure3)],gcomputeN[5,which(names(gcomputeN)==measure3)],gcomputeN[6,which(names(gcomputeN)==measure3)])/bench2)  ##linear delta01
  
  
  bias1=100*(rbind(iptwBias1[1,], gcomputeBias1[1,], aiptwBias1[1,],pencompBias1[1,],
                   iptwBias1[1,], gcomputeBias1[2,],aiptwBias1[2,], pencompBias1[2,],  
                   iptwBias1[2,], gcomputeBias1[1,],aiptwBias1[3,], pencompBias1[3,]))
  
  bias2=100*(rbind(iptwBias2[1,], gcomputeBias2[1,], aiptwBias2[1,],pencompBias2[1,],
                   iptwBias2[1,], gcomputeBias2[2,],aiptwBias2[2,], pencompBias2[2,],  
                   iptwBias2[2,], gcomputeBias2[1,],aiptwBias2[3,], pencompBias2[3,]))
  
  bias1=format(bias1, digit=0)
  bias2=format(bias2, digit=0)
  
  biasAll=cbind(rep("&", dim(bias1)[1]), bias1[,1], rep("&", dim(bias1)[1]), bias1[,2],  rep("&", dim(bias1)[1]), bias1[,3], rep("&", dim(bias1)[1]),
                bias1[,4], rep("&", dim(bias1)[1]), bias1[,5],  rep("&", dim(bias1)[1]), bias1[,6], rep("&", dim(bias1)[1]),
                bias1[,7], rep("&", dim(bias1)[1]), bias1[,8],  rep("&", dim(bias1)[1]), bias1[,9], rep("&", dim(bias1)[1]),rep("&", dim(bias1)[1]),
                bias2[,1], rep("&", dim(bias2)[1]), bias2[,2], rep("&", dim(bias2)[1]), bias2[,3], rep("&", dim(bias2)[1]),
                bias2[,4], rep("&", dim(bias2)[1]), bias2[,5], rep("&", dim(bias2)[1]), bias2[,6], rep("&", dim(bias2)[1]),
                bias2[,7], rep("&", dim(bias2)[1]), bias2[,8], rep("&", dim(bias2)[1]), bias2[,9],rep("\\\\", dim(bias2)[1]))
  
  row.names(biasAll)=c("IPTW(A)","g-computation(A)", "AIPTW(A)", "PENCOMP(A)",
                       "IPTW(A)","g-computation(B)", "AIPTW(B)", "PENCOMP(B)",
                       "IPTW(C)", "g-computation(A)", "AIPTW(C)", "PENCOMP(C)")
  
  write.table(format(biasAll, digits=0), paste(DIREC, "paperTables/sampleSize", sampleSize, "relativeBiasTable.txt",sep=""), sep="", quote=F, col.names=F)
  
  
}
