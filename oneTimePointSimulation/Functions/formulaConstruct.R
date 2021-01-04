


######formula construct
formulaF=function(varList, y.name){
  return ( as.formula(paste(y.name, "~ ", paste(c(varList), collapse = "+"))) )
}


###formula construct for generalized model
formulaGAM=function(varList, y.name, spline){
  return ( as.formula(paste(y.name, " ~ ", paste(c(varList), collapse = "+"), "+", spline)) )
}




