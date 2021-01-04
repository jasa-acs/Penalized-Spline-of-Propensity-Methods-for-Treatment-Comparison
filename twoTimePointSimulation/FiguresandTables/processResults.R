
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
  
  sd11=sd(result[,1]) ###estimated standard error
  width11=mean(abs(result[,4]-result[,3]))  ###mean confidence width
  
  sd11Boot=mean(result[,2])   ###mean bootstrap standard error
  
  finalOut=c( bias11, bias11/truth, sd11, RMSE11, coverage_ps11, width11, dim(result)[1], sd11Boot)  
  names(finalOut)=c( "bias", "biasPercent", "sd", "RMSE", "coverage", "widthCI", "num.sim", "sdBoot")  
  return(finalOut)
  
}




###################for two time point treatment, there are 3 estimates, delta11=E(Y11-Y00), delta10=E(Y10-Y00), and delta01=E(Y01-Y00) ##############
###result-each row has 12 entires, the first 4 entries include estimate, se, and confidence interval for delta11,
### the next 4 entires for delta10, and final 4 entires for delta01
processTotal=function(resultTotal, truth11, truth10, truth01){
  
  
  result11 = resultTotal[ , 1:4]  ###estimate, se and CI for delta11
  result10 = resultTotal[ , (1:4)+(4*1)]  ###delta10
  result01 = resultTotal[ , (1:4)+(4*2)]  ##delta01
  
  out=c(processResult(result = result11, truth = truth11),
        
        processResult(result = result10, truth = truth10),
        
        processResult(result = result01, truth = truth01) )

  return(out)
  
}
