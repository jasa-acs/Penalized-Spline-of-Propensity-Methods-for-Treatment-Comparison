###simulating data sets

#beta1=c(1.5, 1, 0.1)
#beta2=c(1.5, 1, 0.1)
#beta3=c(0.75, 0.5, 0.05)


###for one time point treatment-we consider these factors 1) linear vs nonlinear in covariates;
###2) level of confounding (beta1, beta2, beta3), strength of association between covariates and outcome/treatment
###3) sample size, 200, 500, 1000
###1000 simulations were run for sample size of 500, and 500 simulations for sample size of 200 and 1000
###three covariates L1, L2, L3, outcome Y, binary treatment indicator A0

###this script estimates the truth ATE by generating a large number of times


simulateTruth=function(sampleSize=100000, beta1, beta2, beta3, seed.num, linear="true"){
  
  set.seed(seed.num) 
  
  n <- sampleSize
  
  simdat <- data.frame(L1 = rnorm(n, 0, 1), L2=rnorm(n, 0, 1), L3=rnorm(n, 0, 1)) ###baseline covariate x
  
  ### varying the strength of associaton between x1, x2, and treatment assignment
  a.lin <-beta1*simdat$L1 + beta2*(simdat$L2) + beta3*(simdat$L1 * simdat$L2)
  
  pa <- exp(a.lin)/(1 + exp(a.lin)) ###treatment probability
  #summary(pa)
  
  
  simdat$A0 <- rbinom(n, 1, prob = pa)###treatment status
  mean(simdat$A0)
  
  simdat$L2_sq=simdat$L2*simdat$L2
  simdat$L3_sq=simdat$L3*simdat$L3
  
  if(linear){  ###linear in covariates
    
    simdat$y1 <- 5 + 3*simdat$L2 + simdat$L3 + rnorm(n, 0, 1) ##potential outcome under treatment
    simdat$y0 <- simdat$L2 + simdat$L3 + rnorm(n, 0, 1) ###potential outcome under control
    
  } else {  ###nonlinear in covariates, including quadratic terms
    
    simdat$y1 <- 5 + 3*simdat$L2 + simdat$L3 + 2*simdat$L2_sq + 2*simdat$L3_sq + rnorm(n, 0, 1) ##potential outcome under treatment
    simdat$y0 <- simdat$L2 + simdat$L3 + rnorm(n, 0, 1) ###potential outcome under control
    
  }
  
  return ( mean(simdat$y1-simdat$y0) )  ####true treatment effects by simulating a very large population
  
}

