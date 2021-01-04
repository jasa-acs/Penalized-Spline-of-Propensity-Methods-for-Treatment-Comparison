####this script outputs the models depending on the factors
### In the simulation, we considered these 3 cases-1) both prediction and propensity models correctly specified
### 2) prediction models misspecified, 3) propensity models misspecified
### 
###modelType-correct for case 1, mispred for case 2, and misweight for case 3

modelSpec=function(modelType="correct", linear="true") {
  
  
  ###################################################################################################
  ###propensity models at both time points
  propenVarListA0.Init=NULL
  propenVarListA1.Init=NULL
  
  ###intermediate outcome
  outcomeVarListL1.Init=NULL
  outcomeVarListL1b.Init=NULL
  
  ###final outcome
  outcomeVarListY11.Init=NULL
  outcomeVarListY10.Init=NULL
  outcomeVarListY01.Init=NULL
  outcomeVarListY00.Init=NULL
  
  
  if (linear){  ###if linear in covariates, no interaction terms
    
    if (modelType=="correct"){
      
      propenVarListA0.Init=c("L0","L0b")
      propenVarListA1.Init=c("L1L0", "L1L0:A0", "L1bL0b", "L1bL0b:A0")
      
      outcomeVarListL1.Init=c("L0", "L0b")
      outcomeVarListL1b.Init=c("L1", "L0b")
      
      outcomeVarListY11.Init=c("L0", "L1", "L0b", "L1b")
      outcomeVarListY10.Init=c("L0", "L1", "L0b", "L1b")
      outcomeVarListY01.Init=c("L0", "L1", "L0b", "L1b")
      outcomeVarListY00.Init=c("L0", "L1", "L0b", "L1b")
      
      
    } else if (modelType=="misPred"){
      
      propenVarListA0.Init=c("L0","L0b")
      propenVarListA1.Init=c("L1L0", "L1L0:A0", "L1bL0b", "L1bL0b:A0")
      
      outcomeVarListL1.Init=c("L0", "L0b")
      outcomeVarListL1b.Init=c("L1", "L0b")
      
      outcomeVarListY11.Init=c("L0", "L1")
      outcomeVarListY10.Init=c("L0", "L1")
      outcomeVarListY01.Init=c("L0", "L1")
      outcomeVarListY00.Init=c("L0", "L1")
      
      
    } else if (modelType=="misWeight"){
      
      propenVarListA0.Init=c("L0","L0b")
      propenVarListA1.Init=c("L0", "L1", "L0b")
      
      outcomeVarListL1.Init=c("L0", "L0b")
      outcomeVarListL1b.Init=c("L1", "L0b")
      
      outcomeVarListY11.Init=c("L0", "L1", "L0b", "L1b")
      outcomeVarListY10.Init=c("L0", "L1", "L0b", "L1b")
      outcomeVarListY01.Init=c("L0", "L1", "L0b", "L1b")
      outcomeVarListY00.Init=c("L0", "L1", "L0b", "L1b")
      
    }
    
  } else {  ###nonlinear in covariates
    
    if (modelType=="correct"){
      
      propenVarListA0.Init=c("L0","L0b")
      propenVarListA1.Init=c("L1L0", "L1L0:A0", "L1bL0b", "L1bL0b:A0")
      
      outcomeVarListL1.Init=c("L0", "L0b")
      outcomeVarListL1b.Init=c("L1", "L0b")
      
      outcomeVarListY11.Init=c("L0", "L1", "L0b", "L1b", "L1L1bInt")
      outcomeVarListY10.Init=c("L0", "L1", "L0b", "L1b", "L1L1bInt")
      outcomeVarListY01.Init=c("L0", "L1", "L0b", "L1b", "L1L1bInt")
      outcomeVarListY00.Init=c("L0", "L1", "L0b", "L1b", "L1L1bInt")
      
      
    } else if (modelType=="misPred"){
      
      propenVarListA0.Init=c("L0","L0b")
      propenVarListA1.Init=c("L1L0", "L1L0:A0", "L1bL0b", "L1bL0b:A0")
      
      outcomeVarListL1.Init=c("L0", "L0b")
      outcomeVarListL1b.Init=c("L1", "L0b")
      
      outcomeVarListY11.Init=c("L0", "L1")
      outcomeVarListY10.Init=c("L0", "L1")
      outcomeVarListY01.Init=c("L0", "L1")
      outcomeVarListY00.Init=c("L0", "L1")
      
      
    } else if (modelType=="misWeight"){
      
      propenVarListA0.Init=c("L0","L0b")
      propenVarListA1.Init=c("L0", "L1", "L0b")
      
      outcomeVarListL1.Init=c("L0", "L0b")
      outcomeVarListL1b.Init=c("L1", "L0b")
      
      outcomeVarListY11.Init=c("L0", "L1", "L0b", "L1b", "L1L1bInt")
      outcomeVarListY10.Init=c("L0", "L1", "L0b", "L1b", "L1L1bInt")
      outcomeVarListY01.Init=c("L0", "L1", "L0b", "L1b", "L1L1bInt")
      outcomeVarListY00.Init=c("L0", "L1", "L0b", "L1b", "L1L1bInt")
      
    }
    
  }
  
  return( list(propenVarListA0.Init=propenVarListA0.Init, propenVarListA1.Init=propenVarListA1.Init,
               ###intermediate outcome
               outcomeVarListL1.Init=outcomeVarListL1.Init, outcomeVarListL1b.Init=outcomeVarListL1b.Init,
               
               ###final outcome
               outcomeVarListY11.Init=outcomeVarListY11.Init, outcomeVarListY10.Init=outcomeVarListY10.Init,
               outcomeVarListY01.Init=outcomeVarListY01.Init, outcomeVarListY00.Init=outcomeVarListY00.Init) )
}




