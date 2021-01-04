###simulating a large data set to estimate the truths-theta11, theta10, theta01

#coef1Var=c(-0.5, -0.8, -0.8)
#coef2Var=c(-0.1, -0.1, -0.5)
#coef3Var=c(0.2, 0.6, 1.1)

###for tow-time point treatment-we consider these factors 1) linear vs nonlinear in covariates;
###2) level of confounding (coef1, coef2, coef3), strength of association between covariates and outcome/treatment
###3) sample size, 200, 500, 1000
###1000 simulations were run for sample size of 500, and 500 simulations for sample size of 200 and 1000
###two baseline covariates, L0, L0b, two intermediate covariates L1, L1b, and final outcome Y, binary treatment A0, A1 at both time points
###level-confounding degree: high, moderate, low

simulateTruth=function(sampleSize=100000, level, seed.num, linear){
  
  set.seed(seed.num) 
  
  coef1=coef2=coef3=NULL
  
  if(level=="high"){
    
    coef1=-0.8
    coef2=-0.5
    coef3=1.1
    
  } else if(level=="moderate"){
    
    coef1=-0.8
    coef2=-0.1
    coef3=0.6
    
  } else {
    
    coef1=-0.5
    coef2=-0.1
    coef3=0.2
    
  }
  
  n=sampleSize
  
  L0 = rnorm(n, 0.2, 1)  ###higher the value, more healthy the patient is, so less likely to get treated
  L0b =  rnorm(n, 0.2, 1)  ### side effects, higher the value, more dangerous the drug is for the patient
  
  
  lin1_a =  coef1 * L0 - 0.3*L0b - 0.01   ###0.5, 0.8
  lin1_b = exp(lin1_a) / (1 + exp(lin1_a)) 
  summary(lin1_b)
  A0 = rbinom(n, 1, prob = lin1_b)  ### first time point treatment assignment
  mean(A0)
  
  
  #######################
  L1_0t = rnorm(n, L0 + 0.5*L0b, 1)
  L1_1t = rnorm(n, L0 + 0.5 + 0.5*L0 + 0.5*L0b, 1)
  L1 = L1_0t*(1-A0) + L1_1t*A0
  
  ########################
  L1b_0t = rnorm(n, 0.3*L1_0t + L0b, 1)
  L1b_1t = rnorm(n, 0.4*L1_1t + L0b, 1)
  L1b = L1b_0t*(1-A0) + L1b_1t*A0
  
  cor(L1, L1b)
  
  ########################################
  lin2_a = coef2 * (L1 - L0) + coef3 * (L1 - L0)*A0 + (-0.1) * (L1b - L0b) + coef3 * (L1b - L0b) * A0 - 0.01  ### (L1 - L0) * (1 - A0)
  lin2_b = exp(lin2_a) / (1 + exp(lin2_a))
  summary(lin2_b)
  A1 = rbinom(n, 1, prob = lin2_b)   ###second time point treatment assignment
  mean(A1)
  
  if(linear){
    
    y00 =   5 + 10 * 0 + 10 * 0 + L0 + L1_0t + L0*0 + L1_0t*0 + L0b + L0b * 0 * (0.5)  + L1b_0t + L1b_0t * 0 * (0.5) + rnorm(n, 0, 1)
    y01 =   5 + 10 * 0 + 10 * 1 + L0 + L1_0t + L0*0 + L1_0t*1 + L0b + L0b * 0 * (0.5)  + L1b_0t + L1b_0t * 1 * (0.5) + rnorm(n, 0, 1)
    y10 =   5 + 10 * 1 + 10 * 0 + L0 + L1_1t + L0*1 + L1_1t*0 + L0b + L0b * 1 * (0.5)  + L1b_1t + L1b_1t * 0 * (0.5) + rnorm(n, 0, 1)
    y11 =   5 + 10 * 1 + 10 * 1 + L0 + L1_1t + L0*1 + L1_1t*1 + L0b + L0b * 1 * (0.5)  + L1b_1t + L1b_1t * 1 * (0.5) + rnorm(n, 0, 1)
    
  } else {
    
    y00=    5 + 10 * 0 + 10 * 0 + L0 + L1_0t + L0*0 + L1_0t*0 + L0b + L0b * 0 * (0.5)  + L1b_0t + L1b_0t * 0 * (0.5) + 0.7*L1_0t*L1b_0t + rnorm(n, 0, 1)
    y01 =   5 + 10 * 0 + 10 * 1 + L0 + L1_0t + L0*0 + L1_0t*1 + L0b + L0b * 0 * (0.5)  + L1b_0t + L1b_0t * 1 * (0.5) + 0.8*L1_0t*L1b_0t + rnorm(n, 0, 1)
    y10 =   5 + 10 * 1 + 10 * 0 + L0 + L1_1t + L0*1 + L1_1t*0 + L0b + L0b * 1 * (0.5)  + L1b_1t + L1b_1t * 0 * (0.5) + L1_1t*L1b_1t +  rnorm(n, 0, 1)
    y11 =   5 + 10 * 1 + 10 * 1 + L0 + L1_1t + L0*1 + L1_1t*1 + L0b + L0b * 1 * (0.5)  + L1b_1t + L1b_1t * 1 * (0.5) + 1.6*L1_1t*L1b_1t + rnorm(n, 0, 1)
    
  }
  
  simdat <- data.frame(L0, A0, L1, L1_0t, L1_1t, L1b_0t, L1b_1t, A1, y, y00, y10, y01, y11, L0b, L1b)
  
  return( c( theta11=mean(simdat$y11-simdat$y00), theta10=mean(simdat$y10-simdat$y00), theta01=mean(simdat$y01-simdat$y00) ) )
}

