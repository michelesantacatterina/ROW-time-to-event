########################################################################
########################################################################
#***********************************************************************
#
#   In this file, we provide functions for:
#   1. data generation
#   2. compute hazard ratios
#   3. compute coverage and standard errors
# 
#***********************************************************************
########################################################################
########################################################################


library(CBPS)
library(mvtnorm)
library(cobalt)
library(survival)
library(SuperLearner)
library(gurobi)
library(WeightIt)
library(Matching)
library(sbw)
library(ebal)
# library(optmatch)
library(boot)
library(splines)
library(mfp)
library(survey)



#################    ##################################    #################
#################    ##################################    #################
#
# Data Generation
#
#################    ##################################    #################
#################    ##################################    #################


# baseline hazard: Weibull

# N       = sample size    
# beta    = fixed effect parameter
# rateC   = rate parameter of the exponential distribution of C


#***********************************************************************
#
# This function generates the data for the Simulations section for the 
# binary and continuous treatments;
# I followed Section 5.1 of https://onlinelibrary.wiley.com/doi/full/10.1002/sim.8837
#
#***********************************************************************
data_generation <- function(n,true_haz, 
                            beta_tr, beta_miss, rateC, 
                            type_treatment, scenario)
{
  
  ########################################################################
  #untreated group event time
  T0 <- rexp(n,1)
  
  ########################################################################
  # covariates
  minU <- -.5
  maxU <- 0.5
  Z1 <- Z2 <- runif(n,minU,maxU)
  
  minU <- -1
  maxU <- 1
  Z3 <- Z4 <- runif(n,minU,maxU)
  
  # continuous
  X1 <- -.3 + (0.5*T0)/(T0+1) + 0.4*Z1
  X2 <- -0.3 + log(T0+2) + Z2
  X3 <- 1/(T0+2) + Z1
  X7 <- -.3 + (0.5*T0)/(T0+1) + 0.4*Z3
  X8 <- -0.3 + log(T0+2) + rnorm(n,0,0.5)
  X9 <- 1/(T0+2) + rnorm(n,0,0.5)
  
  # binary
  p4 <- 0.2 + 0.6/(T0+3)
  p5 <- 0.3 + 0.4/(0.5*T0 + 2)
  p6 <- 1/(T0+1)
  X4 <- rbinom(n,1,p4)
  X5 <- rbinom(n,1,p5)
  X6 <- rbinom(n,1,p6)
  
  p10 <- 0.2 + 0.6/(T0+3)
  p11 <- 0.3 + 0.4/(0.5*T0 + 2)
  p12 <- 1/(T0+1)
  X10 <- rbinom(n,1,p10)
  X11 <- rbinom(n,1,p11)
  X12 <- rbinom(n,1,p12)
  
  X <- cbind(X1,X2,X3,X4,X5,X6,X7,X8,X9,X10,X11,X12)

  ########################################################################
  #treatment and time
  
  #************************
  if(type_treatment == "binary"){
    mbeta <- mean(beta_tr * ( 
      -0.1*X[,1] -
        0.3*X[,2] + 
        0.1*X[,3] + 
        0.6*X[,4] + 
        0.4*X[,5] + 
        0.5*X[,6] +
        -0.1*X[,7] - 
        0.3*X[,8] + 
        0.1*X[,9] + 
        0.6*X[,10] + 
        0.4*X[,11] + 
        0.5*X[,12]
    ))
    prt   <- 1./(1. + exp( mbeta - beta_tr * ( 
      -0.1*X[,1] -
        0.3*X[,2] + 
        0.1*X[,3] + 
        0.6*X[,4] + 
        0.4*X[,5] + 
        0.5*X[,6] +
        -0.1*X[,7] - 
        0.3*X[,8] + 
        0.1*X[,9] + 
        0.6*X[,10] + 
        0.4*X[,11] + 
        0.5*X[,12]) 
    )#end exp
    )#end 1/
    
    Tr <- rbinom(n,1,prt) 
    
    Tlat  <- T0 * exp(-log(true_haz)*Tr)
    
  }#end trt binary  
  
  #************************
  if(type_treatment == "continuous"){
    
    meanlog <- -1.5
      if(scenario==1){
        temp <- 0.005
        mTr <- mean( temp*( -0.1*X[,1] - 0.3*X[,2] + 0.1*X[,3] + 0.6*X[,4] + 0.4*X[,5] + 0.5*X[,6] -
                              0.1*X[,7] - 0.3*X[,8] + 0.1*X[,9] + 0.6*X[,10] + 0.4*X[,11] + 0.5*X[,12]) )
        Tr <- - mTr + temp*( -0.1*X[,1] - 0.3*X[,2] + 0.1*X[,3] + 0.6*X[,4] + 0.4*X[,5] + 0.5*X[,6] -
                               0.1*X[,7] - 0.3*X[,8] + 0.1*X[,9] + 0.6*X[,10] + 0.4*X[,11] + 0.5*X[,12]) + rlnorm(n,0,0.1)       #Tr <- - mTr + ( 0.1*X[,1] - 0.3*X[,2] + 0.1*X[,3] + 0.6*X[,4] + 0.4*X[,5] + 0.5*X[,6] )/10 + rnorm(n)   # Weibull latent event times
        
         
        Tlat  <- T0 * exp(-log(true_haz)*Tr)
      }#end scenario 1
    
      if(scenario==2){
        temp <- 0.1
        mTr <- mean( temp*( -0.1*X[,1] - 0.3*X[,2] + 0.1*X[,3] + 0.6*X[,4] + 0.4*X[,5] + 0.5*X[,6] -
                              0.1*X[,7] - 0.3*X[,8] + 0.1*X[,9] + 0.6*X[,10] + 0.4*X[,11] + 0.5*X[,12]) )
        Tr <- - mTr + temp*( -0.1*X[,1] - 0.3*X[,2] + 0.1*X[,3] + 0.6*X[,4] + 0.4*X[,5] + 0.5*X[,6] -
                               0.1*X[,7] - 0.3*X[,8] + 0.1*X[,9] + 0.6*X[,10] + 0.4*X[,11] + 0.5*X[,12]) + rlnorm(n,0,0.1)       #Tr <- - mTr + ( 0.1*X[,1] - 0.3*X[,2] + 0.1*X[,3] + 0.6*X[,4] + 0.4*X[,5] + 0.5*X[,6] )/10 + rnorm(n)   # Weibull latent event times
        
        Tlat  <- T0 * exp(-log(true_haz)*Tr)
      }#end scenario 2
      
      if(scenario==3){
        temp <- 0.05
        mTr <- mean( temp*( -0.1*X[,1] - 0.3*X[,2] + 0.1*X[,3] + 0.6*X[,4] + 0.4*X[,5] + 0.5*X[,6] -
                              0.1*X[,7] - 0.3*X[,8] + 0.1*X[,9] + 0.6*X[,10] + 0.4*X[,11] + 0.5*X[,12]) )
        Tr <- - mTr + temp*( -0.1*X[,1] - 0.3*X[,2] + 0.1*X[,3] + 0.6*X[,4] + 0.4*X[,5] + 0.5*X[,6] -
                               0.1*X[,7] - 0.3*X[,8] + 0.1*X[,9] + 0.6*X[,10] + 0.4*X[,11] + 0.5*X[,12]) + rlnorm(n,0,0.1)       #Tr <- - mTr + ( 0.1*X[,1] - 0.3*X[,2] + 0.1*X[,3] + 0.6*X[,4] + 0.4*X[,5] + 0.5*X[,6] )/10 + rnorm(n)   # Weibull latent event times
        
        
        Tlat  <- T0 * exp(-log(true_haz)*Tr)
      }#end scenario 3
      
      
  }#end trt continuous
  
  
  
  ########################################################################
  #censoring times
  C <- rexp(n=n, rate=rateC)
  
  
  
  ########################################################################
  # follow-up times and event indicators
  time      <- pmin(Tlat, C)
  status    <- as.numeric(Tlat <= C)
  censoring <- 1-status
  
  ########################################################################
  # misspecification
  
  if(type_treatment == "binary"){
    W1  <- X2/(1+exp(X1))+1
    W2  <- (X7*X8/25+1)^3
    W3 <-  (X2+X3+2)^2
    W7 <- exp(X7/2)
    W8 <- (X7*X8/25+1)^3
    W9 <- (X8+X9)^2

    X1  <- X1*(1-beta_miss)+W1*(beta_miss)
    X2  <- X2*(1-beta_miss)+W2*(beta_miss)
    X3  <- X3*(1-beta_miss)+W3*(beta_miss)
  }
  
  if(type_treatment == "continuous"){

    W1  <- X2/(1+exp(X1))+1
    W2  <- (X7*X8/25+1)^3
    W3 <-  (X2+X3+2)^2
    W7 <- exp(X7/2)
    W8 <- (X7*X8/25+1)^3
    W9 <- (X8+X9)^2

    X1  <- X1*(1-beta_miss)+W1*(beta_miss)
    X2  <- X2*(1-beta_miss)+W2*(beta_miss)
    X3  <- X3*(1-beta_miss)+W3*(beta_miss)

  }
  
  
  X <- cbind(X1,X2,X3,X4,X5,X6,X7,X8,X9,X10,X11,X12)
  ########################################################################
  # data set
  dta <- data.frame(id=1:n,
                    time=time,
                    status=status,
                    censoring=censoring,
                    Tr=Tr,
                    X=X)
  dta <- dta[order(dta$time,-dta$status),]
  return(dta)
  
} # end data_generation





#***********************************************************************
#
# This function generates the data for the Simulations section in the supplementary material for the 
# binary and continuous treatments when increasing the sample size and number of covariates
#
#***********************************************************************


data_generation_cova <- function(n, lambda, rho, delta, 
                            beta_tr, beta_miss, beta_ou, rateC, 
                            type_treatment, scenario, num_cova)
{

  Sigma <- diag(num_cova)
  X     <- mvrnorm(n , rep(0, num_cova), Sigma)
  
  ########################################################################
  #treatment and time
  
  #************************
  if(type_treatment == "binary"){
    mbeta <- mean( rep(beta_tr,num_cova)%*%t(X) )
    prt   <- 1./(1. + exp( mbeta - t(rep(beta_tr,num_cova)%*%t(X)) ))
    
    Tr <- rbinom(n,1,prt) 
    
    v     <- runif(n=n)
    Tlat  <- (- log(v) / (lambda * exp(Tr*delta + t(rep(beta_ou,num_cova)%*%t(X)) )))^(1 / rho)
    
  }#end trt binary  
  
  #************************
  if(type_treatment == "continuous"){
    if(scenario==1){
      mTr <- mean( rep(beta_tr,num_cova)%*%t(X) )
      Tr <- - mTr + (t(rep(beta_tr,num_cova)%*%t(X))) + rnorm(n)   # Weibull latent event times
      v     <- runif(n=n)
      Tlat  <- (- log(v) / (lambda * exp(Tr*delta + t(rep(beta_ou,num_cova)%*%t(X)) )))^(1 / rho)
    }#end scenario 1
    
    if(scenario==2){
      temp <- 0.5
      mTr <- mean( rep(beta_tr,num_cova)%*%t(X) )
      Tr <- - mTr + (t(rep(beta_tr,num_cova)%*%t(X))) + rlnorm(n,0,temp)
      v     <- runif(n=n)
      Tlat  <- (- log(v) / (lambda * exp(Tr*delta + t(rep(beta_ou,num_cova)%*%t(X)) )))^(1 / rho)
    }#end scenario 2
    
  }#end trt continuous
  
  
  ########################################################################
  #censoring times
  C <- rexp(n=n, rate=rateC)
  
  
  
  ########################################################################
  # follow-up times and event indicators
  time      <- pmin(Tlat, C)
  status    <- as.numeric(Tlat <= C)
  censoring <- 1-status
  
  # data set
  dta <- data.frame(id=1:n,
                    time=time,
                    status=status,
                    censoring=censoring,
                    Tr=Tr,
                    X=X)
  dta <- dta[order(dta$time,-dta$status),]
  return(dta)
  
} # end data_generationNum cov





#################    ##################################    #################
#################    ##################################    #################
#
# ROW
#
#################    ##################################    #################
#################    ##################################    #################

#***********************************************************************
#
# This function provides ROW weights using Gurobi
#
#***********************************************************************

row <- function(confounders,intervention,delta=0.01){
  
  Xs <- scale(cbind(confounders))
  trs <- as.numeric(scale(intervention))
  
  tol <- 1e-08
  reptimes <- dim(Xs)[2]
  model <- list()
  model$A <- matrix(c(rep(Xs*trs,2),rep(1,n)), nrow = (reptimes*2+1), byrow = T)
  model$rhs <- c(rep(delta,reptimes),rep(-delta,reptimes),1)
  model$modelsense <- "min"
  model$Q <- diag(n)
  model$obj <- rep(1/(n-1),n)
  model$sense <- c(rep("<=",reptimes),rep(">=",reptimes),"=")
  model$lb <- rep(tol, n)
  model$vtypes <- "C"
  params <- list(Presolve = 2, OutputFlag = 0)
  res <- gurobi(model,params)
  return(res)
  
}





#################    ##################################    #################
#################    ##################################    #################
#
# IPW WEIGHTS
#
#################    ##################################    #################
#################    ##################################    #################


#***********************************************************************
#
# This function provides a logistic model with splines 
#
#***********************************************************************

SL.glm.BS <- function(Y, X, newX, family, obsWeights, model = TRUE, ...) 
{
  if (is.matrix(X)) {
    X = as.data.frame(X)
  }

  
  m1 <- median(X[,1])
  m2 <- median(X[,2])
  m3 <- median(X[,3])
  m7 <- median(X[,7])
  m8 <- median(X[,8])
  m9 <- median(X[,9])
  
  
  fit.glm <- glm(Y ~ bs(X.X1,knots=m1) + 
                     bs(X.X2,knots=m2) + 
                     bs(X.X3,knots=m3) +
                     bs(X.X7,knots=m7) + 
                     bs(X.X8,knots=m8) + 
                     bs(X.X9,knots=m9) 
                     , data = X, family = family, weights = obsWeights, 
                 model = model)
  if (is.matrix(newX)) {
    newX = as.data.frame(newX)
  }
  pred <- predict(fit.glm, newdata = newX, type = "response")
  fit <- list(object = fit.glm)
  class(fit) <- "SL.glm.BS"
  out <- list(pred = pred, fit = fit)
  return(out)
}



SL.glm.FP <- function(Y, X, newX, family, obsWeights, model = TRUE, ...) 
{
  if (is.matrix(X)) {
    X = as.data.frame(X)
  }
  
  
  
  fit.glm <- glm(Y ~ fp(X.X1,df=2) + 
                   fp(X.X2,df=2) + 
                   fp(X.X3,df=2) +
                   fp(X.X7,df=2) + 
                   fp(X.X8,df=2) + 
                   fp(X.X9,df=2) 
                 , data = X, family = family, weights = obsWeights, 
                 model = model)
  if (is.matrix(newX)) {
    newX = as.data.frame(newX)
  }
  pred <- predict(fit.glm, newdata = newX, type = "response")
  fit <- list(object = fit.glm)
  class(fit) <- "SL.glm.FP"
  out <- list(pred = pred, fit = fit)
  return(out)
}

#***********************************************************************
#
# This function provides IPW weights by using Superlearner
#
#***********************************************************************
gps_weights <- function(dta, type_treatment, sl.lib_gps, sl.lib_ps, stabilize, truncate, qt){
  ##################
  # PS/GPS Estimation
  ##################
  
  if(type_treatment == "binary"){
    
    #Balancing covariates between treatment groups (binary)
    # temp  <- dim(dta)[2]
    # Y     <- dta$Tr
    # X     <- dta[,(temp-5):temp]
    # SL    <- SuperLearner(Y, X ,SL.library = sl.lib_ps, family = "binomial")
    # ps    <- SL$SL.predict
    # ipw   <- ( dta$Tr/ps + (1-dta$Tr)/(1-ps)  )
    ps <- tryCatch( weightit(Tr ~ X.X1 +X.X2+X.X3 +X.X4 +X.X5 +X.X6 +
                               X.X7 +X.X8 +X.X9 +X.X10 +X.X11 +X.X12, data = dta,
                    method = "super", stabilize = stabilize,
                    SL.library = sl.lib_ps, family="binomial")
                    , error=function(e) NULL)#endtryCatch
    
    ipw <- ps$weights
    
    if(truncate==T & is.null(ps)!=T){
      indicator_trunc <- tryCatch( (ps$weights>quantile(ps$weights,qt))*(ps$weights<quantile(ps$weights,1-qt)) , error=function(e) NULL)#endtryCatch
      ipw <- indicator_trunc*ps$weights
      ipw[which(ipw==0)] <- 1e-08
    }#end if truncate
    
    return(list(ps_w = ipw))
    
  }
  else{
    
    gps <- tryCatch( weightit(Tr ~ X.X1 +X.X2+X.X3 +X.X4 +X.X5 +X.X6 +
                                X.X7 +X.X8 +X.X9 +X.X10 +X.X11 +X.X12, data = dta,
                   method = "super", density = "dnorm",
                   stabilize = stabilize,
                   SL.library = sl.lib_gps, family="gaussian")
                  , error=function(e) NULL)#endtryCatch
    
    ipw <- gps$weights
    
    if(truncate==T & is.null(gps)!=T){
      indicator_trunc <- (gps$weights>quantile(gps$weights,qt))*(gps$weights<quantile(gps$weights,1-qt))
      ipw <- indicator_trunc*gps$weights
      ipw[which(ipw==0)] <- 1e-08
    }#end if truncate
  
    
    return(list(ps_w = ipw))
    
  }
  

}#end if(method=="IPW")




#################    ##################################    #################
#################    ##################################    #################
#
# COMPUTE HAZARD RATIO
#
#################    ##################################    #################
#################    ##################################    #################

#***********************************************************************
#
# This function compute the hazard ratio for all methods defined in 'method'
#
#***********************************************************************

compute_haz_hat <- function(itera,method,rateCens,type_treatment,beta_pos,beta_miss,n,scenario,delta){
  
  errors <- rep(0,itera)
  coef  <- mse <- time <- mcens <- timet <- NA
  bal   <- matrix(rep(NA,itera*12),ncol=12,nrow=itera)
  
  for(i in 1:itera){
    #print(i)
    if(i==1 | i%%25==0){print(paste("Iteration",i,"out of",itera))}  
    dta <- read.csv(paste("/Users/santam13/Documents/NYU/ROW/ROW-time-to-event/code/simulations/data_tmp/",
                          rateCens,"_type_treatment_",type_treatment,"_lackover_",beta_pos,"_misspe_",beta_miss,"_scenario_",
                          scenario,"_n_",n,"_",i,".csv",sep=""))
    
    confounders   <- data.frame(dta$X.X1,dta$X.X2,dta$X.X3,dta$X.X4,dta$X.X5,dta$X.X6,
                                dta$X.X7,dta$X.X8,dta$X.X9,dta$X.X10,dta$X.X11,dta$X.X12)
    intervention  <- dta$Tr
    
    A <- matrix(scale(cbind(confounders))*as.numeric(scale(intervention)), nrow = dim(confounders)[2], byrow = T)
    
    ####################################################################################################
    if(method=="ROW"){
      timet     <-  tryCatch( system.time( row_w <- row(confounders,intervention,delta)$x)
                    , error=function(e) NULL)#endtryCatch
      
      if(is.null(timet)!=TRUE){
        fit     <- tryCatch( summary(coxph(Surv(time, status) ~ Tr, data=dta, weights=row_w))
                    , error=function(e) NULL)#endtryCatch
        
        if(is.null(fit)!=TRUE){
          coef[i] <- fit$coefficients[1]
          mse[i]  <- (coef[i]-log(true_haz))^2
          if(type_treatment=="binary"){
            bal[i,] <- abs(bal.tab(intervention~confounders,
                                   weights=row_w,
                                   method="weighting",
                                   s.d.denom="pooled")$Balance$Diff.Adj)
          }#end if(type_treatment=="binary")
          else{
            # bal[i,] <- abs(bal.tab(intervention~confounders,
            #                        weights=row_w,
            #                        method="weighting")$Balance$Corr.Adj)
            
            bal[i,] <- abs(A%*%row_w)
            
          }#end else 
        }#end if NULL
        
      }#end if NULL
      else{
        errors[i] <- 1
      }#end else
        
    }#end row
    
    
    ####################################################################################################
    if(method=="CBPS"){
      timet    <-  tryCatch( system.time( cbps_w  <- CBPS(scale(Tr) ~  scale(X.X1) + scale(X.X2) + scale(X.X3) + 
                                                                       scale(X.X4) + scale(X.X5) + scale(X.X6) + 
                                                                       scale(X.X7) + scale(X.X8) + scale(X.X9) + 
                                                                       scale(X.X10) + scale(X.X11) + scale(X.X12),
                                                          method = "exact", data=dta, ATT=0)$weights)
                            , error=function(e) NULL)#endtryCatch
      
      if(is.null(timet)!=TRUE){
        fit     <- tryCatch( summary(coxph(Surv(time, status) ~ Tr, data=dta, weights=cbps_w))
                             , error=function(e) NULL)#endtryCatch
        
        if(is.null(fit)!=TRUE){
          coef[i] <- fit$coefficients[1]
          mse[i]  <- (coef[i]-log(true_haz))^2
          
          if(type_treatment=="binary"){
            bal[i,] <- abs(bal.tab(intervention~confounders,
                                   weights=cbps_w,
                                   method="weighting",
                                   s.d.denom="pooled")$Balance$Diff.Adj)
          }#end if(type_treatment=="binary")
          else{
            # bal[i,] <- abs(bal.tab(intervention~confounders,
            #                        weights=cbps_w,
            #                        method="weighting")$Balance$Corr.Adj)
            
            bal[i,] <- abs(A%*%cbps_w)
          }#end else 
          
        }#end if NULL
      }#end if NULL
      else{
        errors[i] <- 1
      }#end else
    }#end cbps
    
    
    ####################################################################################################
    if(method=="npCBPS"){
      timet    <-  tryCatch( system.time( capture.output(npcbps_w  <- npCBPS(scale(Tr) ~ scale(X.X1) + scale(X.X2) + scale(X.X3) + 
                                                                               scale(X.X4) + scale(X.X5) + scale(X.X6) + 
                                                                               scale(X.X7) + scale(X.X8) + scale(X.X9) + 
                                                                               scale(X.X10) + scale(X.X11) + scale(X.X12),
                                                                             data=dta,print.level=0)$weights))
                            , error=function(e) NULL)#endtryCatch
      
      if(is.null(timet)!=TRUE){
        fit     <- tryCatch( summary(coxph(Surv(time, status) ~ Tr, data=dta, weights=npcbps_w))
                             , error=function(e) NULL)#endtryCatch
        
        if(is.null(fit)!=TRUE){
          coef[i] <- fit$coefficients[1]
          mse[i]  <- (coef[i]-log(true_haz))^2
          
          if(type_treatment=="binary"){
            bal[i,] <- abs(bal.tab(intervention~confounders,
                                   weights=npcbps_w,
                                   method="weighting",
                                   s.d.denom="pooled")$Balance$Diff.Adj)
          }#end if(type_treatment=="binary")
          else{
            # bal[i,] <- abs(bal.tab(intervention~confounders,
            #                        weights=npcbps_w,
            #                        method="weighting")$Balance$Corr.Adj)
            
            bal[i,] <- abs(A%*%npcbps_w)
          }#end else 
          
        }#end if NULL
      }#end if NULL
      else{
        errors[i] <- 1
      }#end else
    }#end npCBPS
    
    
    ####################################################################################################
    if(method=="IPW-SL"){
      
      stabilize <- FALSE
      truncate <- FALSE
      qt <- 0.1
      
      timet    <- tryCatch( system.time(gps_Slw <- gps_weights(dta, type_treatment, sl.lib_gps, sl.lib_ps, stabilize, truncate, qt))
                            , error=function(e) NULL)#endtryCatch
      
      if(is.null(timet)!=TRUE){
        fit     <- tryCatch( summary(coxph(Surv(time, status) ~ Tr, data=dta, weights=gps_Slw$ps_w))
                             , error=function(e) NULL)#endtryCatch
        
        if(is.null(fit)!=TRUE){
          coef[i] <- fit$coefficients[1]
          mse[i]  <- (coef[i]-log(true_haz))^2
          
          if(type_treatment=="binary"){
            bal[i,] <- abs(bal.tab(intervention~confounders,
                                   weights=gps_Slw$ps_w,
                                   method="weighting",
                                   s.d.denom="pooled")$Balance$Diff.Adj)
          }#end if(type_treatment=="binary")
          else{
            # bal[i,] <- abs(bal.tab(intervention~confounders,
            #                        weights=gps_Slw$ps_w,
            #                        method="weighting")$Balance$Corr.Adj)
            
            testBal <- tryCatch( abs(A%*%(gps_Slw$ps_w/sum(gps_Slw$ps_w))), error=function(e) NULL)
            if(is.null(testBal)!=TRUE){
              bal[i,] <- testBal
            }
          }#end else 
          
        }#end if NULL
      }#end if NULL
      else{
        errors[i] <- 1
      }#end else
    }#end IPW-SL
    
    
    
    ####################################################################################################
    if(method=="IPW-SL-stab"){
      
      stabilize <- TRUE
      truncate <- FALSE
      qt <- 0.1
      
      timet    <- tryCatch( system.time(gps_Slw <- gps_weights(dta, type_treatment, sl.lib_gps, sl.lib_ps, stabilize, truncate, qt))
                            , error=function(e) NULL)#endtryCatch
      
      if(is.null(timet)!=TRUE){
        fit     <- tryCatch( summary(coxph(Surv(time, status) ~ Tr, data=dta, weights=gps_Slw$ps_w))
                             , error=function(e) NULL)#endtryCatch
        
        if(is.null(fit)!=TRUE){
          coef[i] <- fit$coefficients[1]
          mse[i]  <- (coef[i]-log(true_haz))^2
          
          if(type_treatment=="binary"){
            bal[i,] <- abs(bal.tab(intervention~confounders,
                                   weights=gps_Slw$ps_w,
                                   method="weighting",
                                   s.d.denom="pooled")$Balance$Diff.Adj)
          }#end if(type_treatment=="binary")
          else{
            # bal[i,] <- abs(bal.tab(intervention~confounders,
            #                        weights=gps_Slw$ps_w,
            #                        method="weighting")$Balance$Corr.Adj)
            
            testBal <- tryCatch( abs(A%*%(gps_Slw$ps_w/sum(gps_Slw$ps_w))), error=function(e) NULL)
            if(is.null(testBal)!=TRUE){
              bal[i,] <- testBal
            }
          }#end else 
          
        }#end if NULL
      }#end if NULL
      else{
        errors[i] <- 1
      }#end else
    }#end IPW-SL-stab
    
    
    ####################################################################################################
    if(method=="IPW-SL-trunc"){
      
      stabilize <- FALSE
      truncate <- TRUE
      qt <- 0.01
      
      timet    <- tryCatch( system.time(gps_Slw <- gps_weights(dta, type_treatment, sl.lib_gps, sl.lib_ps, stabilize, truncate, qt))
                            , error=function(e) NULL)#endtryCatch
      
      if(is.null(timet)!=TRUE){
        fit     <- tryCatch( summary(coxph(Surv(time, status) ~ Tr, data=dta, weights=gps_Slw$ps_w))
                             , error=function(e) NULL)#endtryCatch
        
        if(is.null(fit)!=TRUE){
          coef[i] <- fit$coefficients[1]
          mse[i]  <- (coef[i]-log(true_haz))^2
          
          if(type_treatment=="binary"){
            testBal <- tryCatch( abs(bal.tab(intervention~confounders,
                                   weights=gps_Slw$ps_w,
                                   method="weighting",
                                   s.d.denom="pooled")$Balance$Diff.Adj), error=function(e) NULL)
            if(is.null(testBal)!=TRUE){
              bal[i,] <- testBal
            }
          }#end if(type_treatment=="binary")
          else{
            # bal[i,] <- abs(bal.tab(intervention~confounders,
            #                        weights=gps_Slw$ps_w,
            #                        method="weighting")$Balance$Corr.Adj)
            
            testBal <- tryCatch( abs(A%*%(gps_Slw$ps_w/sum(gps_Slw$ps_w))), error=function(e) NULL)
            if(is.null(testBal)!=TRUE){
              bal[i,] <- testBal
            }
          }#end else 
          
        }#end if NULL
      }#end if NULL
      else{
        errors[i] <- 1
      }#end else
    }#end IPW-SL-stab
    
    
    ####################################################################################################
    if(method=="IPW-BS"){
      
      m1 <- median(dta$X.X1)
      m2 <- median(dta$X.X2)
      m3 <- median(dta$X.X3)

      m7 <- median(dta$X.X7)
      m8 <- median(dta$X.X8)
      m9 <- median(dta$X.X9)
      
      if(type_treatment=="binary"){

        timet    <- tryCatch( system.time(fit_ipwBS <- glm(Tr ~ bs(X.X1,knots=m1) + 
                                                             bs(X.X2,knots=m2) + 
                                                             bs(X.X3,knots=m3) + 
                                                             X.X4 + 
                                                             X.X5 + 
                                                             X.X6 +
                                                             bs(X.X7,knots=m7) + 
                                                             bs(X.X8,knots=m8) + 
                                                             bs(X.X9,knots=m9) + 
                                                             X.X10 + 
                                                             X.X11 + 
                                                             X.X12
                                                             ,
                                                              data = dta,
                                                             family = "binomial") )
                              , error=function(e) NULL)#endtryCatch
        
        if(is.null(timet)!=TRUE){
          pr1 <- fit_ipwBS$fitted.values
          ipwBS_weights <- ( dta$Tr/pr1 + (1-dta$Tr)/(1-pr1)  )
        }
        
      }#end if treatment binary
      
      if(type_treatment=="continuous"){
        
        timet    <- tryCatch( system.time(fit_ipwBS <- weightit(Tr ~ bs(X.X1,knots=m1) + 
                                                                  bs(X.X2,knots=m2) + 
                                                                  bs(X.X3,knots=m3) + 
                                                                  X.X4 + 
                                                                  X.X5 + 
                                                                  X.X6 +
                                                                  bs(X.X7,knots=m7) + 
                                                                  bs(X.X8,knots=m8) + 
                                                                  bs(X.X9,knots=m9) + 
                                                                  X.X10 + 
                                                                  X.X11 + 
                                                                  X.X12,
                                                             data = dta,
                                                             method = "ps",
                                                             density = "dnorm") )
                              , error=function(e) NULL)#endtryCatch
        
        if(is.null(timet)!=TRUE){
          ipwBS_weights <- fit_ipwBS$weights
        }
        
      }#end else
      
      
      if(is.null(timet)!=TRUE){
        
        
        fit     <- tryCatch( summary(coxph(Surv(time, status) ~ Tr, data=dta, weights=ipwBS_weights))
                             , error=function(e) NULL)#endtryCatch
        
        if(is.null(fit)!=TRUE){
          coef[i] <- fit$coefficients[1]
          mse[i]  <- (coef[i]-log(true_haz))^2
          
          if(type_treatment=="binary"){
            bal[i,] <- abs(bal.tab(intervention~confounders,
                                   weights=ipwBS_weights,
                                   method="weighting",
                                   s.d.denom="pooled")$Balance$Diff.Adj)
          }#end if(type_treatment=="binary")
          else{
            # bal[i,] <- abs(bal.tab(intervention~confounders,
            #                        weights=ipwBS_weights,
            #                        method="weighting")$Balance$Corr.Adj)
            
            testBal <- tryCatch( abs(A%*%(ipwBS_weights/sum(ipwBS_weights))), error=function(e) NULL)
            if(is.null(testBal)!=TRUE){
              bal[i,] <- testBal
            }
          }#end else 
          
        }#end if NULL
      }#end if NULL
      else{
        errors[i] <- 1
      }#end else
    }#end IPW-Splines
    
    
    
    ####################################################################################################
    if(method=="IPW-FP"){
      
      if(type_treatment=="binary"){
        
        timet    <- tryCatch( system.time(fit_ipwBS <- glm(Tr ~ fp(X.X1,df=2) + 
                                                             fp(X.X2,df=2) + 
                                                             fp(X.X3,df=2) + 
                                                             X.X4 + 
                                                             X.X5 + 
                                                             X.X6 +
                                                             fp(X.X7,df=2) + 
                                                             fp(X.X8,df=2) + 
                                                             fp(X.X9,df=2) + 
                                                             X.X10 + 
                                                             X.X11 + 
                                                             X.X12
                                                           ,
                                                           data = dta,
                                                           family = "binomial") )
                              , error=function(e) NULL)#endtryCatch
        
        if(is.null(timet)!=TRUE){
          pr1 <- fit_ipwBS$fitted.values
          ipwBS_weights <- ( dta$Tr/pr1 + (1-dta$Tr)/(1-pr1)  )
        }
        
      }#end if treatment binary
      
      if(type_treatment=="continuous"){
        
        timet    <- tryCatch( system.time(fit_ipwBS <- weightit(Tr ~ fp(X.X1,df=2) + 
                                                                  fp(X.X2,df=2) + 
                                                                  fp(X.X3,df=2) + 
                                                                  X.X4 + 
                                                                  X.X5 + 
                                                                  X.X6 +
                                                                  fp(X.X7,df=2) + 
                                                                  fp(X.X8,df=2) + 
                                                                  fp(X.X9,df=2) + 
                                                                  X.X10 + 
                                                                  X.X11 + 
                                                                  X.X12,
                                                                data = dta,
                                                                method = "ps",
                                                                density = "dnorm") )
                              , error=function(e) NULL)#endtryCatch
        
        if(is.null(timet)!=TRUE){
          ipwBS_weights <- fit_ipwBS$weights
        }
        
      }#end else
      
      
      if(is.null(timet)!=TRUE){
        
        
        fit     <- tryCatch( summary(coxph(Surv(time, status) ~ Tr, data=dta, weights=ipwBS_weights))
                             , error=function(e) NULL)#endtryCatch
        
        if(is.null(fit)!=TRUE){
          coef[i] <- fit$coefficients[1]
          mse[i]  <- (coef[i]-log(true_haz))^2
          
          if(type_treatment=="binary"){
            bal[i,] <- abs(bal.tab(intervention~confounders,
                                   weights=ipwBS_weights,
                                   method="weighting",
                                   s.d.denom="pooled")$Balance$Diff.Adj)
          }#end if(type_treatment=="binary")
          else{
            # bal[i,] <- abs(bal.tab(intervention~confounders,
            #                        weights=ipwBS_weights,
            #                        method="weighting")$Balance$Corr.Adj)
            
            testBal <- tryCatch( abs(A%*%(ipwBS_weights/sum(ipwBS_weights))), error=function(e) NULL)
            if(is.null(testBal)!=TRUE){
              bal[i,] <- testBal
            }
          }#end else 
          
        }#end if NULL
      }#end if NULL
      else{
        errors[i] <- 1
      }#end else
    }#end IPW-Splines
    
    ####################################################################################################
    if(method=="BalSL"){
      if(type_treatment=="binary"){
        SL.library <- sl.lib_ps
        stop.method <- "es.mean"
        
        
        timet    <- tryCatch( system.time(ps_BalSL <- weightit(Tr ~ X.X1 + X.X2 + X.X3 + X.X4 + X.X5 + X.X6 +
                                                                 X.X7 + X.X8 + X.X9 + X.X10 + X.X11 + X.X12,
                                                               data = dta,
                                                               method = "super",
                                                               estimand = 'ATE',
                                                               SL.library =  SL.library ,
                                                               SL.method = "method.balance",
                                                               stop.method = stop.method) )
                              , error=function(e) NULL)#endtryCatch
        
      }#end if treatment binary
      
      if(type_treatment=="continuous"){
        SL.library <- sl.lib_gps
        stop.method <- "p.mean"
        
        timet    <- tryCatch( system.time(ps_BalSL <- weightit(Tr ~ X.X1 + X.X2 + X.X3 + X.X4 + X.X5 + X.X6 +
                                                                 X.X7 + X.X8 + X.X9 + X.X10 + X.X11 + X.X12,
                                                               data = dta,
                                                               method = "super",
                                                               density = "dnorm",
                                                               SL.library =  SL.library ,
                                                               SL.method = "method.balance",
                                                               stop.method = stop.method) )
                              , error=function(e) NULL)#endtryCatch
        
        
      }#end else
    
      
      if(is.null(timet)!=TRUE){
        fit     <- tryCatch( summary(coxph(Surv(time, status) ~ Tr, data=dta, weights=ps_BalSL$weights))
                             , error=function(e) NULL)#endtryCatch
        
        if(is.null(fit)!=TRUE){
          coef[i] <- fit$coefficients[1]
          mse[i]  <- (coef[i]-log(true_haz))^2
          
          if(type_treatment=="binary"){
            bal[i,] <- abs(bal.tab(intervention~confounders,
                                   weights=ps_BalSL$weights,
                                   method="weighting",
                                   s.d.denom="pooled")$Balance$Diff.Adj)
          }#end if(type_treatment=="binary")
          else{
            # bal[i,] <- abs(bal.tab(intervention~confounders,
            #                        weights=ps_BalSL$weights,
            #                        method="weighting")$Balance$Corr.Adj)
            
            testBal <- tryCatch( abs(A%*%(ps_BalSL$weights/sum(ps_BalSL$weights))), error=function(e) NULL)
            if(is.null(testBal)!=TRUE){
              bal[i,] <- testBal
            }
          }#end else 
          
        }#end if NULL
      }#end if NULL
      else{
        errors[i] <- 1
      }#end else
    }#end Balance Super Learner 
    
    
    
    ####################################################################################################
    if(method=="GBM"){
      if(type_treatment=="binary"){
        stop.method <- "es.mean"
        
       
        timet    <- tryCatch( system.time(ps_gbm <- weightit(Tr ~ X.X1 + X.X2 + X.X3 + X.X4 + X.X5 + X.X6 +
                                                               X.X7 + X.X8 + X.X9 + X.X10 + X.X11 + X.X12,
                                                               data = dta,
                                                               method = "gbm",
                                                               estimand = 'ATE',
                                                               stop.method = stop.method) )
                              , error=function(e) NULL)#endtryCatch
        
      }#end if treatment binary
      
      if(type_treatment=="continuous"){
        stop.method <- "p.mean"
        
        timet    <- tryCatch( system.time(ps_gbm <- weightit(Tr ~ X.X1 + X.X2 + X.X3 + X.X4 + X.X5 + X.X6 +
                                                               X.X7 + X.X8 + X.X9 + X.X10 + X.X11 + X.X12,
                                                               data = dta,
                                                               method = "gbm",
                                                               density = "dnorm",
                                                               stop.method = stop.method) )
                              , error=function(e) NULL)#endtryCatch
        
        
      }#end else
      
      
      if(is.null(timet)!=TRUE){
        fit     <- tryCatch( summary(coxph(Surv(time, status) ~ Tr, data=dta, weights=ps_gbm$weights))
                             , error=function(e) NULL)#endtryCatch
        
        if(is.null(fit)!=TRUE){
          coef[i] <- fit$coefficients[1]
          mse[i]  <- (coef[i]-log(true_haz))^2
          
          if(type_treatment=="binary"){
            bal[i,] <- abs(bal.tab(intervention~confounders,
                                   weights=ps_gbm$weights,
                                   method="weighting",
                                   s.d.denom="pooled")$Balance$Diff.Adj)
          }#end if(type_treatment=="binary")
          else{
            # bal[i,] <- abs(bal.tab(intervention~confounders,
            #                        weights=ps_gbm$weights,
            #                        method="weighting")$Balance$Corr.Adj)
            
            testBal <- tryCatch( abs(A%*%(ps_gbm$weights/sum(ps_gbm$weights))), error=function(e) NULL)
            if(is.null(testBal)!=TRUE){
              bal[i,] <- testBal
            }
          }#end else 
          
        }#end if NULL
      }#end if NULL
      else{
        errors[i] <- 1
      }#end else
    }#end GBM
    
    
    
    ####################################################################################################
    if(method=="SBW"){
      #SBW
      balsbw         <-  list()
      balsbw$bal_cov <- c("X.X1","X.X2","X.X3","X.X4","X.X5","X.X6",
                          "X.X7","X.X8","X.X9","X.X10","X.X11","X.X12")
      balsbw$bal_tol <- 0.01
      balsbw$bal_std <- "group"
      
      timet <- tryCatch( system.time({capture.output( out0 <-  sbw(dta, ind = "Tr",  bal = balsbw, 
                                                   sol = list(sol_nam = "gurobi"), par = list(par_est = "ate"))) 
                                      sbw <- out0$dat_weights$sbw_weights
                          })#end system.time 
      , error=function(e) NULL)#endtryCatch
      
      fit <- tryCatch(  coxph(Surv(time, status) ~ Tr , data=dta, weights=sbw)  
                 , error=function(e) NULL)#endtryCatch
      
      if(is.null(timet)!=TRUE){
          if(is.null(fit)!=TRUE){
            coef[i] <- fit$coefficients[1]
            mse[i]  <- (coef[i]-log(true_haz))^2
            #print(i)
            bal[i,] <- abs(bal.tab(intervention~confounders,
                                 weights=sbw,
                                 method="weighting",
                                 s.d.denom="pooled")$Balance$Diff.Adj)
        }#end if null fit
      }#end if null
      else{
        errors[i] <- 1
      }#end else
      
    }#end method == "SBW"
    
    
    ####################################################################################################
    if(method == "PSM"){
      
      stabilize <- FALSE
      truncate <- FALSE
      qt <- 0.1
      
      #compute ipw weight
      gps_Slw <- tryCatch( gps_weights(dta, type_treatment, sl.lib_gps, sl.lib_ps, stabilize, truncate, qt)
                            , error=function(e) NULL)#endtryCatch
      
      timet <- tryCatch(  system.time({ m1 <- Match(Tr = intervention, X = gps_Slw$ps_w, estimand = "ATE")
      })#end system.time
      , error=function(e) NULL)#endtryCatch
      
      if(is.null(timet)!=TRUE){
        
        matches <- factor(rep(m1$index.treated, 2))
        matchedsample <- cbind(matches, dta[c(m1$index.control, m1$index.treated),])
        fit <- tryCatch(  coxph(Surv(time, status) ~ Tr + strata(matches), data=matchedsample)
                   , error=function(e) NULL)#endtryCatch
        
        # fit <- tryCatch(  coxph(Surv(time, status) ~ Tr, data=mdata)  
        #             , error=function(e) NULL)#endtryCatch
        
        if(is.null(fit)!=TRUE){
          coef[i] <- fit$coefficients[1]
          mse[i]  <- (coef[i]-log(true_haz))^2
          bal[i,] <- abs(bal.tab(m1,intervention~confounders)$Balance$Diff.Adj)
          #bal[i,] <- abs(bal.tab(m2, disp.subclass = TRUE)$Balance$Diff.Adj[-1])
        }#end if null fit
        
      }#end if null
      else{
        errors[i] <- 1
      }#end else

      
    }#end if method == PSM
    
    
    ####################################################################################################
    if(method == "EBAL"){
      
      #compute ipw weight
      timet <- tryCatch(  system.time({ebal   <- weightit(Tr ~ X.X1 +X.X2 +X.X3 +X.X4 +X.X5 +X.X6  +
                                                            X.X7 + X.X8 + X.X9 + X.X10 + X.X11 + X.X12, 
                                                          data = dta,
                                                          method = "ebal", 
                                                          estimand = "ATE")
      })#end system.time
      , error=function(e) NULL)#endtryCatch
      
      if(is.null(timet)!=TRUE){
        fit <- tryCatch(  coxph(Surv(time, status) ~ Tr , data=dta, weights=ebal$w)  
                        , error=function(e) NULL)#endtryCatch
      
        if(is.null(fit)!=TRUE){
        coef[i] <- fit$coefficients[1]
        mse[i]  <- (coef[i]-log(true_haz))^2
        bal[i,] <- abs(bal.tab(intervention~confounders,
                               weights=ebal$w,
                               method="weighting",
                               s.d.denom="pooled")$Balance$Diff.Adj)
        }#end if null
      }#end if null
      else{
        errors[i] <- 1
      }#end else
    
    }#end if method == EBAL
    
    
    ####################################################################################################
    if(method=="OM"){
      timet        <- tryCatch( system.time(fit       <- summary(coxph(Surv(time, status) ~ Tr + X.X1 + X.X2 + X.X3 + X.X4 + X.X5 + X.X6  +
                                                                         X.X7 + X.X8 + X.X9 + X.X10 + X.X11 + X.X12
                                                                       , data=dta)))
                                , error=function(e) NULL)#endtryCatch

      if(is.null(timet)!=TRUE){
        coef[i]   <- fit$coefficients[1]
        mse[i]    <- (coef[i]-true_haz)^2
        bal[i,]   <- rep(0,6)
      }#end if NULL
      else{
        errors[i] <- 1
      }#end else
    }#end OM
    
    
    ####################################################################################################
    if(method=="naive"){
      #naive
      timet      <- tryCatch( system.time(fit     <- summary(coxph(Surv(time, status) ~ Tr, data=dta)))
                              , error=function(e) NULL)#endtryCatch

      if(is.null(timet)!=TRUE){
        coef[i] <- fit$coefficients[1]
        mse[i]  <- (coef[i]-log(true_haz))^2
        
        if(type_treatment=="binary"){
          bal[i,] <- abs(bal.tab(intervention~confounders, s.d.denom = 'pooled')$Balance$Diff.Un)
        }#end if(type_treatment=="binary")
        else{
          bal[i,] <- abs(bal.tab(intervention~confounders)$Balance$Corr.Un)
        }#end else 
        

      }#end if NULL
      else{
        errors[i] <- 1
      }#end else
    }#end method=="naive
    
    if(is.null(timet)!=T){
      time[i]   <- timet[1] + timet[2]
    }
    mcens[i]  <- mean(dta$censoring,na.rm=T)
    
  }#end itera
  
  Mcens   <- mean(mcens,na.rm=T)
  
  bias    <- mean(abs(coef-log(true_haz)),na.rm=T)
  Mse     <- sqrt(mean(mse,na.rm=T))
  Mbal    <- apply(bal,2,mean,na.rm=T)
  Mtime   <- mean(time,na_rm=T)
  
  if(is.null(bias)!=T & is.null(Mse)!=T & is.null(Mbal)!=T & is.null(Mtime)!=T ){
    return(list(bias=bias,mse=Mse,bal=Mbal,time=Mtime,pcens=Mcens,errors=errors))
  }
  else{return(list(bias=NA,mse=NA,bal=NA,time=NA))}
  
}#end compute_haz_hat




compute_haz_hat_cova <- function(itera,method,rateCens,type_treatment,beta_pos,beta_miss,n,scenario,delta,num_cova){
  
  coef  <- mse <- time <- mcens <- timet <- NA
  bal   <- matrix(rep(NA,itera*num_cova),ncol=num_cova,nrow=itera)
  
  for(i in 1:itera){
    #print(i)
    if(i==1 | i%%25==0){print(paste("Iteration",i,"out of",itera))}  
    dta <- read.csv(paste("/Users/santam13/Documents/NYU/ROW/ROW-time-to-event/code/simulations/data_tmp/",
                          rateCens,"_type_treatment_",type_treatment,"_lackover_",beta_pos,"_misspe_",beta_miss,"_scenario_",
                          scenario,"_n_",n,"_","num_cova",num_cova,"_",i,".csv",sep=""))
    
    confounders   <- dta[,6:dim(dta)[2]]
    intervention  <- dta$Tr
    
    ####################################################################################################
    if(method=="ROW"){
      timet     <-  tryCatch( system.time( row_w <- row(confounders,intervention,delta)$x)
                              , error=function(e) NULL)#endtryCatch
      
      if(is.null(timet)!=TRUE){
        fit     <- tryCatch( summary(coxph(Surv(time, status) ~ Tr, data=dta, weights=row_w))
                             , error=function(e) NULL)#endtryCatch
        
        if(is.null(fit)!=TRUE){
          coef[i] <- fit$coefficients[1]
          mse[i]  <- (coef[i]-log(true_haz))^2
          if(type_treatment=="binary"){
            bal[i,] <- tryCatch( abs(bal.tab(intervention~confounders,
                                   weights=row_w,
                                   method="weighting",
                                   s.d.denom="pooled")$Balance$Diff.Adj), error=function(e) NULL)#endtryCatch
          }#end if(type_treatment=="binary")
          else{
            bal[i,] <- tryCatch( abs(bal.tab(intervention~confounders,
                                   weights=row_w,
                                   method="weighting")$Balance$Corr.Adj), error=function(e) NULL)#endtryCatch
          }#end else 
        }#end if NULL
        
      }#end if NULL
    }#end row
    
    if(is.null(timet)!=T){
      time[i]   <- timet[1] + timet[2]
    }
    mcens[i]  <- mean(dta$censoring,na.rm=T)
    
  }#end itera
  
  Mcens   <- mean(mcens,na.rm=T)
  
  bias    <- (mean(coef,na.rm=T)-true_haz)
  Mse     <- mean(mse,na.rm=T)
  Mbal    <- apply(bal,2,mean,na.rm=T)
  Mtime   <- mean(time,na_rm=T)
  
  if(is.null(bias)!=T & is.null(Mse)!=T & is.null(Mbal)!=T & is.null(Mtime)!=T ){
    return(list(bias=bias,mse=Mse,bal=Mbal,time=Mtime,pcens=Mcens))
  }
  else{return(list(bias=NA,mse=NA,bal=NA,time=NA))}

}#end compute_haz_hat_cova




#################    ##################################    #################
#################    ##################################    #################
#
# FUNCTIONS FOR COVERAGE
#
#################    ##################################    #################
#################    ##################################    #################


#***********************************************************************
#
# This function provides bootstrap standard errors
#
#***********************************************************************
boot_row <- function(data, indices, cor.type){
  haz_hatb <-NA
  resampled_sample  <- data[indices,]
  confounders_s       <- data.frame(resampled_sample$X.X1,resampled_sample$X.X2,
                                  resampled_sample$X.X3,resampled_sample$X.X4,
                                  resampled_sample$X.X5,resampled_sample$X.X6,
                                  resampled_sample$X.X7,resampled_sample$X.X8,
                                  resampled_sample$X.X9,resampled_sample$X.X10,
                                  resampled_sample$X.X11,resampled_sample$X.X12)
  intervention_s  <- resampled_sample$Tr
  
  #ROW
  row_w <- tryCatch( row_w <- row(confounders_s,intervention_s)$x,error=function(e) NULL)
  
  if(is.null(row_w)!=TRUE){
    #estimate hazard ratio
    fitb      <- coxph(Surv(time, status) ~ Tr, data=resampled_sample, weights=row_w)
    if(is.null(fitb)!=TRUE){
      haz_hatb  <- fitb$coefficients[1]
    }#end if NULL
  }#end if(res NULL)
  exp(haz_hatb)
}#end boot_row



bootW_row <- function(data, indices, cor.type){
  haz_hatb <-NA
  resampled_sample  <- data[indices,]

    #estimate hazard ratio
    fitb      <- coxph(Surv(time, status) ~ Tr, data=resampled_sample)
    if(is.null(fitb)!=TRUE){
      haz_hatb  <- fitb$coefficients[1]
  }#end if(res NULL)
  exp(haz_hatb)
}#end boot_row



boot_ipw <- function(data, indices, cor.type){
  haz_hatb <-NA
  resampled_sample  <- data[indices,]

  stabilize <- FALSE
  truncate <- FALSE
  qt <- 0.1
  
  #IPW
  
  sl.lib_gps2 <- c('SL.glm')
  sl.lib_ps2  <- c('SL.glm')
  
  timet    <- tryCatch( system.time(gps_Slw <- gps_Slw <- gps_weights(resampled_sample, type_treatment, sl.lib_gps2, sl.lib_ps2, stabilize, truncate, qt))
                        , error=function(e) NULL)#endtryCatch
  
  if(is.null(timet)!=TRUE){
    fit     <- tryCatch( summary(coxph(Surv(time, status) ~ Tr, data=resampled_sample, weights=gps_Slw$ps_w))
                         , error=function(e) NULL)#endtryCatch
    if(is.null(fit)!=TRUE){
      haz_hatb  <- fit$coefficients[1]
    }#end if NULL
  }#end if(res NULL)
  exp(haz_hatb)
}#end boot_ipw


bootW_ipw <- function(data, indices, cor.type){
  haz_hatb <-NA
  resampled_sample  <- data[indices,]
  
  #estimate hazard ratio
  fitb      <- coxph(Surv(time, status) ~ Tr, data=resampled_sample)
  if(is.null(fitb)!=TRUE){
    haz_hatb  <- fitb$coefficients[1]
  }#end if(res NULL)
  exp(haz_hatb)
}#end boot_row


#***********************************************************************
#
# This function provides standard errors as described in the Practical guidelines section
#
#***********************************************************************
compute_errors <- function(itera,method,rateCens,type_treatment,beta_pos,beta_miss,n,scenario,itera_boot,true_haz,delta){
  
  sen <- ser <- seb <- coefn <- sebw <- NA
  timeR <- timeN <- timeB <- timeBW <- NA
  cov95r <- cov95b <- cov95n <- cov95bw <-  cov95b_n <- cov95bw_n <- NA
  mcens <- Mcens <- NA
  
  for(i in 1:itera){
    #print(i)
    if(i==1 | i%%25==0){print(paste("Iteration",i,"out of",itera))}  
    #print(i)
    dta <- read.csv(paste("/Users/santam13/Documents/NYU/ROW/ROW-time-to-event/code/simulations/data_tmp/",
                          rateCens,"_type_treatment_",type_treatment,"_lackover_",beta_pos,"_misspe_",beta_miss,"_scenario_",
                          scenario,"_n_",n,"_",i,".csv",sep=""))
    
    confounders   <- data.frame(dta$X.X1,dta$X.X2,dta$X.X3,dta$X.X4,dta$X.X5,dta$X.X6,
                                dta$X.X7,dta$X.X8,dta$X.X9,dta$X.X10,dta$X.X11,dta$X.X12)
    intervention  <- dta$Tr
  
    ####################################################################################################
    if(method=="ROW"){
      
      row_w <- tryCatch( row(confounders,intervention,delta)$x, error=function(e) NULL)#endtryCatch
      
      if(is.null(row_w)!=TRUE){
          ############
          #robust
          surveyDesign  <- svydesign(ids=1:n,weights=~row_w,data=dta)
          timer <- system.time({
            fitr <- tryCatch( svycoxph(Surv(time, status) ~ Tr,design=surveyDesign) 
                                                  , error=function(e) NULL)#endtryCatch
          
            if(is.null(fitr)!=TRUE){
              ser[i]    <- sqrt(fitr$var)
              capture <- capture.output( tempS     <- summary(fitr) )
              cov95r[i] <- (tempS$conf.int[3]<=true_haz)*(true_haz<=tempS$conf.int[4])
            }#end if NULL
          })#end system timer
          
          ############
          #naive
          timen <- system.time({
            #I need to make weights sum to N when using robust=F
            fitn          <- tryCatch( coxph(Surv(time, status) ~ Tr, data=dta, weights=row_w*n, robust=F) , error=function(e) NULL)#endtryCatch
          
            if(is.null(fitn)!=TRUE){
              coefn[i]  <- fitn$coefficients[1]
              sen[i]    <- sqrt(fitn$var)
              tempS     <- summary(fitn)
              cov95n[i] <- (tempS$conf.int[3]<=true_haz)*(true_haz<=tempS$conf.int[4])
            }#end if NULL
          })#end system timer
          
      }#end if NULL
      
      ############
      #bootstrap
      timeb <- system.time({
        myBootstrap <- boot(dta, boot_row, R=itera_boot)
        seb[i]      <- sd(myBootstrap$t)
        ci          <- boot.ci(myBootstrap,conf=0.95,type=c("perc","norm"))
        cov95b[i]   <- (ci$percent[4]<=true_haz)*(true_haz<=ci$percent[5])
        cov95b_n[i] <- (ci$normal[2]<=true_haz)*(true_haz<=ci$normal[3])
      })#end system timer
      
      
      ############
      #Weighted bootstrap
      timebw <- system.time({
        myBootstrapW <- boot(dta, bootW_row, R=itera_boot, weights=row_w)
        sebw[i]      <- sd(myBootstrapW$t)
        ci           <- boot.ci(myBootstrapW,conf=0.95,type=c("perc","norm"))
        cov95bw[i]   <- (ci$percent[4]<=true_haz)*(true_haz<=ci$percent[5])
        cov95bw_n[i] <- (ci$normal[2]<=true_haz)*(true_haz<=ci$normal[3])
      })#end system timer
      
      timeR[i] <- timer[1] +  timer[2]
      timeN[i] <- timen[1] +  timen[2]
      timeB[i] <- timeb[1] +  timeb[2]
      timeBW[i] <- timebw[1] +  timebw[2]
      
      mcens[i]  <- mean(dta$censoring,na.rm=T)
      Mcens   <- mean(mcens,na.rm=T)
      
      
    }#end method row
    
    
    ####################################################################################################
    if(method=="IPW"){
      
      
      sl.lib_gps2 <- c('SL.glm')
      sl.lib_ps2  <- c('SL.glm')
      
      gps_Slw <- gps_weights(dta, type_treatment, sl.lib_gps2, sl.lib_ps2, stabilize, truncate, qt)
      
      if(is.null(gps_Slw)!=TRUE){
        ############
        #robust
        surveyDesign  <- svydesign(ids=1:n,weights=~gps_Slw$ps_w,data=dta)
        timer <- system.time({
          fitr <- tryCatch( svycoxph(Surv(time, status) ~ Tr,design=surveyDesign) 
                            , error=function(e) NULL)#endtryCatch
          
          if(is.null(fitr)!=TRUE){
            ser[i]    <- sqrt(fitr$var)
            capture <- capture.output( tempS     <- summary(fitr) )
            cov95r[i] <- (tempS$conf.int[3]<=true_haz)*(true_haz<=tempS$conf.int[4])
          }#end if NULL
        })#end system timer
        
        ############
        #naive
        timen <- system.time({
          fitn          <- tryCatch( coxph(Surv(time, status) ~ Tr, data=dta, weights=gps_Slw$ps_w, robust=F) , error=function(e) NULL)#endtryCatch
          
          if(is.null(fitn)!=TRUE){
            coefn[i]  <- fitn$coefficients[1]
            sen[i]    <- sqrt(fitn$var)
            tempS     <- summary(fitn)
            cov95n[i] <- (tempS$conf.int[3]<=true_haz)*(true_haz<=tempS$conf.int[4])
          }#end if NULL
        })#end system timer
        
      }#end if NULL
      
      ############
      #bootstrap
      timeb <- system.time({
        myBootstrap <- boot(dta, boot_ipw, R=itera_boot)
        seb[i]      <- sd(myBootstrap$t)
        ci          <- boot.ci(myBootstrap,conf=0.95,type=c("perc","norm"))
        cov95b[i]   <- (ci$percent[4]<=true_haz)*(true_haz<=ci$percent[5])
        cov95b_n[i] <- (ci$normal[2]<=true_haz)*(true_haz<=ci$normal[3])
      })#end system timer
      
      
      ############
      #Weighted bootstrap
      timebw <- system.time({
        myBootstrapW <- boot(dta, bootW_ipw, R=itera_boot, weights=gps_Slw$ps_w)
        sebw[i]      <- sd(myBootstrapW$t)
        ci           <- boot.ci(myBootstrapW,conf=0.95,type=c("perc","norm"))
        cov95bw[i]   <- (ci$percent[4]<=true_haz)*(true_haz<=ci$percent[5])
        cov95bw_n[i] <- (ci$normal[2]<=true_haz)*(true_haz<=ci$normal[3])
      })#end system timer
      
      timeR[i] <- timer[1] +  timer[2]
      timeN[i] <- timen[1] +  timen[2]
      timeB[i] <- timeb[1] +  timeb[2]
      timeBW[i] <- timebw[1] +  timebw[2]
      
      mcens[i]  <- mean(dta$censoring,na.rm=T)
      Mcens   <- mean(mcens,na.rm=T)
      
      
    }#end method row
    
    
  }#end for (itera)
  
  if(is.null(coefn)!=T & is.null(sen)!=T & is.null(ser)!=T & is.null(seb)!=T  ){
    return(list(ses=sd(coefn,na.rm=T),
                sen=mean(sen,na.rm=T),
                seb=mean(seb,na.rm=T),
                ser=mean(ser,na.rm=T),
                sebw=mean(sebw,na.rm=T),
                timeR=mean(timeR,na.rm=T),
                timeN=mean(timeN,na.rm=T),
                timeB=mean(timeB,na.rm=T),
                timeBW=mean(timeBW,na.rm=T),
                cov95r=mean(cov95r,na.rm=T),
                cov95b=mean(cov95b,na.rm=T),
                cov95n=mean(cov95n,na.rm=T),
                cov95b_n=mean(cov95b_n,na.rm=T),
                cov95bw=mean(cov95bw,na.rm=T),
                cov95bw_n=mean(cov95bw_n,na.rm=T),
                pcens=Mcens
                )
           )
  }
  else{return(list(coefn=NA,
                   sen=NA,
                   ser=NA,
                   seb=NA,
                   timeR=NA,
                   timeN=NA,
                   timeB=NA,
                   cov95b=NA,
                   cov95r=NA,
                   cov95n=NA,
                   pcens=NA))}
  
}#end function compute errors








