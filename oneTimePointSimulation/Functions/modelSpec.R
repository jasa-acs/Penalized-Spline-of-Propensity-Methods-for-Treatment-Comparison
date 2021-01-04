####this script outputs the models depending on the factors
### In the simulation, we considered these 3 cases-1) both prediction and propensity models correctly specified
### 2) prediction models misspecified, 3) propensity models misspecified
### 
###modelType-correct for case 1, mispred for case 2, and misweight for case 3

modelSpec=function(modelType="correct", linear="true") {
  
  
  ###################################################################################################
  propenVarList.Init=NULL  ###propensity model
  outcomeVarList0.Init=NULL ##model for Y0
  outcomeVarList1.Init=NULL  ##model for Y1
  
  if (linear){  ###if linear in covariates
    
    if (modelType=="correct"){
      
      propenVarList.Init=c("L1", "L2", "L1L2")
      outcomeVarList0.Init=c("L2", "L3")
      outcomeVarList1.Init=c("L2", "L3")
      
    } else if (modelType=="misPred"){
      
      propenVarList.Init=c("L1", "L2", "L1L2")
      outcomeVarList0.Init=c("L3")
      outcomeVarList1.Init=c("L3")
      
    } else if (modelType=="misWeight"){
      
      propenVarList.Init=c("L1")
      outcomeVarList0.Init=c("L2", "L3")
      outcomeVarList1.Init=c("L2", "L3")
      
    }
    
  } else {  ###nonlinear in covariates
    
    if (modelType=="correct"){
      
      propenVarList.Init=c("L1", "L2", "L1L2")
      outcomeVarList0.Init=c("L2", "L3")
      outcomeVarList1.Init=c("L2", "L3", "L2_sq", "L3_sq")
      
    } else if (modelType=="misPred"){
      
      propenVarList.Init=c("L1", "L2", "L1L2")
      outcomeVarList0.Init=c("L3")
      outcomeVarList1.Init=c("L3")
      
    } else if (modelType=="misWeight"){
      
      propenVarList.Init=c("L1")
      outcomeVarList0.Init=c("L2", "L3")
      outcomeVarList1.Init=c("L2", "L3", "L2_sq", "L3_sq")
      
    }
    
  }
  
  
  return( list(propen.model=propenVarList.Init, modely0=outcomeVarList0.Init, modely1=outcomeVarList1.Init) )
}

