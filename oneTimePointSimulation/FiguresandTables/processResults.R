
############################################################################################################
###result is a matrix, with each row consisting of estimate, sd, and confidence interval
###result-results from simAllMethod.R
###each row in result has 4 entries, estimate, se, and 95% confidence interval lower and upper
###
processResult=function(result, truth) {
  
  result=result[!is.na(result[,1]),]
  
  ###coverage rate
  total_ps11=numeric(dim(result)[1])
  for (g in 1:dim(result)[1]) {
    total_ps11[g]=as.numeric(result[g,3] <= truth & result[g,4] >= truth)
  }
  coverage_ps11=sum(total_ps11)/length(total_ps11)
  coverage_ps11
  
  ###bias and RMSE
  bias11=mean(result[,1]-(truth))
  estimate11=mean(result[,1])
  temp11=(result[,1]-(truth))^2
  RMSE11=sqrt(mean(temp11))
  
  sd11=sd(result[,1])
  width11=mean(abs(result[,4]-result[,3]))
  
  sd11Boot=mean(result[,2]) 
  
  finalOut=c(bias11, bias11/truth, sd11, RMSE11, coverage_ps11, width11, dim(result)[1], sd11Boot)  
  names(finalOut)=c( "bias", "biasPercent", "sd", "RMSE", "coverage", "widthCI", "num.sim", "sdBoot")  
  return(finalOut)
  
}






