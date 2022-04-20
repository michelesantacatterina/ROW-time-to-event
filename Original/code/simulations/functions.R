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
library(optmatch)
library(boot)



#################    ##################################    #################
#################    ##################################    #################
#
# Data Generation
#
#################    ##################################    #################
#################    ##################################    #################


# baseline hazard: Weibull

# N       = sample size    
# lambda  = scale parameter in h0()
# rho     = shape parameter in h0()
# beta    = fixed effect parameter
# rateC   = rate parameter of the exponential distribution of C


#***********************************************************************
#
# This function generates the data for the Simulations section for the 
# binary and continuous treatments;
#
#***********************************************************************
data_generation <- function(n, lambda, rho, delta, 
                            beta_tr, beta_miss, beta_ou, rateC, 
                            type_treatment, scenario)
{
  
  ########################################################################
  # covariates
  X1 <- rnorm(n,.1,1)
  X2 <- rnorm(n,.1,1)
  X3 <- rlnorm(n,0,.5)
  X4 <- 5*rbeta(n,3,1)
  
  X5 <- sample( 1:4, n, replace=TRUE, prob=c(0.35, 0.25, 0.05, 0.35) )
  X6 <- sample( 0:1, n, replace=TRUE, prob=c(0.75, 0.25) )
  
  X <- cbind(X1,X2,X3,X4,X5,X6)
  

  ########################################################################
  #treatment and time
  
  #************************
  if(type_treatment == "binary"){
    mbeta <- mean(beta_tr*(X[,1] + X[,2] + X[,3] + X[,4] + X[,5]) + X[,6])
    prt   <- 1./(1. + exp( mbeta - beta_tr * (X[,1] + X[,2] + X[,3] + X[,4] + X[,5]) + X[,6]))
    
    Tr <- rbinom(n,1,prt) 
    
    v     <- runif(n=n)
    Tlat  <- (- log(v) / (lambda * exp(Tr*delta + X[,2] + beta_ou*(X[,4] + X[,5] ) + X[,6] )))^(1 / rho)
    
  }#end trt binary  
  
  #************************
  if(type_treatment == "continuous"){
      if(scenario==1){
        temp <- 0.1
        mTr <- mean((X[,1] + X[,2] + X[,3] + X[,4] + X[,5] + X[,6]) )
        Tr <- - mTr + (X[,1] + X[,2] + X[,3] + X[,4] + X[,5] + X[,6]) + rnorm(n)   # Weibull latent event times
        v     <- runif(n=n)
        Tlat  <- (- log(v) / (lambda * exp(Tr*delta + X[,2] + beta_ou*(X[,4] + X[,5] ) + X[,6] )))^(1 / rho)
      }#end scenario 1
    
      if(scenario==2){
        temp <- 0.5
        mTr <- mean((X[,1] + X[,2] + X[,3] + X[,4] + X[,5] + X[,6]) )
        Tr <- - mTr + (X[,1] + X[,2] + X[,3] + X[,4] + X[,5] + X[,6]) + rlnorm(n,0,temp)
        v     <- runif(n=n)
        Tlat  <- (- log(v) / (lambda * exp(Tr*delta + X[,2] + beta_ou*(X[,4] + X[,5] ) + X[,6] )))^(1 / rho)
      }#end scenario 2
      
      if(scenario==3){
        temp <-0.6
        mTr <- mean((X[,1] + X[,2] + X[,3] + X[,4] + X[,5] + X[,6]))
        Tr <- - mTr + (X[,1] + X[,2] + X[,3] + X[,4] + X[,5] + X[,6]) + rlnorm(n,0,temp)
        v     <- runif(n=n)
        Tlat  <- (- log(v) / (lambda * exp(Tr*delta + X[,2] + beta_ou*(X[,4] + X[,5] ) + X[,6] )))^(1 / rho)
      }#end scenario 2
      
      if(scenario==4){
        temp <- 0.7
        mTr <- mean(temp*(X[,1] + X[,2] + X[,3] + X[,4] + X[,5] + X[,6]) )
        Tr <- - mTr + temp*(X[,1] + X[,2] + X[,3] + X[,4] + X[,5] + X[,6]) + rlnorm(n,0,temp)
        v     <- runif(n=n)
        Tlat  <- (- log(v) / (lambda * exp(Tr*delta + X[,2] + beta_ou*(X[,4] + X[,5] ) + X[,6] )))^(1 / rho)
      }#end scenario 4
    
      if(scenario==5){
        temp <- 0.9
        mTr <- mean((X[,1] + X[,2] + X[,3] + X[,4] + X[,5] + X[,6]) )
        Tr <- - mTr + (X[,1] + X[,2] + X[,3] + X[,4] + X[,5] + X[,6]) + rlnorm(n,0,temp)
        v     <- runif(n=n)
        Tlat  <- (- log(v) / (lambda * exp(Tr*delta + X[,2] + beta_ou*(X[,4] + X[,5] ) + X[,6] )))^(1 / rho)
      }#end scenario 5
      
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
    W1  <- (X1+0.5)^2
    W2 <- ((X1*X2)/5 + 1)^2
    W3  <- exp(X3/2)
    W4  <- X4*(1+exp(X3)) + 1
    X1  <- X1*(1-beta_miss)+W1*(beta_miss)
    X2  <- X2*(1-beta_miss)+W2*(beta_miss)
    X3  <- X3*(1-beta_miss)+W3*(beta_miss)
    X4  <- X4*(1-beta_miss)+W4*(beta_miss)
  }
  
  if(type_treatment == "continuous"){
    W1  <- exp(X1/2)
    W2  <- X2*(1+exp(X1)) + 1
    W3 <-  (X1*X3/25 + 0.2)^3
    W4 <-  2*log(abs(X4))
    X1  <- X1*(1-beta_miss)+W1*(beta_miss)
    X2  <- X2*(1-beta_miss)+W2*(beta_miss)
    X3  <- X3*(1-beta_miss)+W3*(beta_miss)
    X4  <- X4*(1-beta_miss)+W4*(beta_miss)
  }
  
  
  X   <- cbind(X1,X2,X3,X4,X5,X6)
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
# This function provides IPW weights by using Superlearner
#
#***********************************************************************
gps_weights <- function(dta, type_treatment, sl.lib_gps, sl.lib_ps){
  ##################
  # PS/GPS Estimation
  ##################
  
  if(type_treatment == "binary"){
    
    #Balancing covariates between treatment groups (binary)
    temp  <- dim(dta)[2]
    Y     <- dta$Tr
    X     <- dta[,(temp-5):temp]
    SL    <- SuperLearner(Y, X ,SL.library = sl.lib_ps, family = "binomial")
    ps    <- SL$SL.predict
    ipw   <- ( dta$Tr/ps + (1-dta$Tr)/(1-ps)  )
    
    return(list(ps = ps, gps_w = ipw))
    
  }
  else{
    
    gps <- weightit(Tr ~ X.X1 +X.X2 +X.X3 +X.X4 +X.X5 +X.X6 , data = dta,
                   method = "super", density = "dnorm",
                   SL.library = sl.lib_gps, family="gaussian")
    #print(max(gps$weights))
    return(list(gps_w = gps$weights))
    
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
  
  coef  <- mse <- time <- mcens <- timet <- NA
  bal   <- matrix(rep(NA,itera*6),ncol=6,nrow=itera)
  
  for(i in 1:itera){
    #print(i)
    if(i==1 | i%%25==0){print(paste("Iteration",i,"out of",itera))}  
    dta <- read.csv(paste("~/path",
                          rateCens,"_type_treatment_",type_treatment,"_lackover_",beta_pos,"_misspe_",beta_miss,"_scenario_",
                          scenario,"_n_",n,"_",i,".csv",sep=""))
    
    confounders   <- data.frame(dta$X.X1,dta$X.X2,dta$X.X3,dta$X.X4,dta$X.X5,dta$X.X6)
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
          mse[i]  <- (coef[i]-true_haz)^2
          if(type_treatment=="binary"){
            bal[i,] <- abs(bal.tab(intervention~confounders,
                                   weights=row_w,
                                   method="weighting",
                                   s.d.denom="pooled")$Balance$Diff.Adj)
          }#end if(type_treatment=="binary")
          else{
            bal[i,] <- abs(bal.tab(intervention~confounders,
                                   weights=row_w,
                                   method="weighting")$Balance$Corr.Adj)
          }#end else 
        }#end if NULL
        
      }#end if NULL
    }#end row
    
    
    ####################################################################################################
    if(method=="CBPS"){
      timet    <-  tryCatch( system.time( cbps_w  <- CBPS(scale(Tr) ~ scale(X.X1) + scale(X.X2) + scale(X.X3) + 
                                                 scale(X.X4) + scale(X.X5) + scale(X.X6),method = "exact", data=dta, ATT=0)$weights)
                            , error=function(e) NULL)#endtryCatch
      
      if(is.null(timet)!=TRUE){
        fit     <- tryCatch( summary(coxph(Surv(time, status) ~ Tr, data=dta, weights=cbps_w))
                             , error=function(e) NULL)#endtryCatch
        
        if(is.null(fit)!=TRUE){
          coef[i] <- fit$coefficients[1]
          mse[i]  <- (coef[i]-true_haz)^2
          
          if(type_treatment=="binary"){
            bal[i,] <- abs(bal.tab(intervention~confounders,
                                   weights=cbps_w,
                                   method="weighting",
                                   s.d.denom="pooled")$Balance$Diff.Adj)
          }#end if(type_treatment=="binary")
          else{
            bal[i,] <- abs(bal.tab(intervention~confounders,
                                   weights=cbps_w,
                                   method="weighting")$Balance$Corr.Adj)
          }#end else 
          
        }#end if NULL
      }#end if NULL
    }#end cbps
    
    
    ####################################################################################################
    if(method=="npCBPS"){
      timet    <-  tryCatch( system.time( capture.output(npcbps_w  <- npCBPS(scale(Tr) ~ scale(X.X1) + scale(X.X2) + scale(X.X3) + 
                                                                    scale(X.X4) + scale(X.X5) + scale(X.X6),data=dta,print.level=0)$weights))
                            , error=function(e) NULL)#endtryCatch
      
      if(is.null(timet)!=TRUE){
        fit     <- tryCatch( summary(coxph(Surv(time, status) ~ Tr, data=dta, weights=npcbps_w))
                             , error=function(e) NULL)#endtryCatch
        
        if(is.null(fit)!=TRUE){
          coef[i] <- fit$coefficients[1]
          mse[i]  <- (coef[i]-true_haz)^2
          
          if(type_treatment=="binary"){
            bal[i,] <- abs(bal.tab(intervention~confounders,
                                   weights=npcbps_w,
                                   method="weighting",
                                   s.d.denom="pooled")$Balance$Diff.Adj)
          }#end if(type_treatment=="binary")
          else{
            bal[i,] <- abs(bal.tab(intervention~confounders,
                                   weights=npcbps_w,
                                   method="weighting")$Balance$Corr.Adj)
          }#end else 
          
        }#end if NULL
      }#end if NULL
    }#end npCBPS
    
    
    ####################################################################################################
    if(method=="IPW"){
      timet    <- tryCatch( system.time(gps_Slw <- gps_weights(dta, type_treatment, sl.lib_gps, sl.lib_ps))
                            , error=function(e) NULL)#endtryCatch
      
      if(is.null(timet)!=TRUE){
        if(type_treatment=="binary"){
        #I want to save the propensity score for PS matching so I don't have to compute it twice
          ps[i,]      <- gps_Slw$ps
        }
        fit     <- tryCatch( summary(coxph(Surv(time, status) ~ Tr, data=dta, weights=gps_Slw$gps_w))
                             , error=function(e) NULL)#endtryCatch
        
        if(is.null(fit)!=TRUE){
          coef[i] <- fit$coefficients[1]
          mse[i]  <- (coef[i]-true_haz)^2
          
          if(type_treatment=="binary"){
            bal[i,] <- abs(bal.tab(intervention~confounders,
                                   weights=gps_Slw$gps_w,
                                   method="weighting",
                                   s.d.denom="pooled")$Balance$Diff.Adj)
          }#end if(type_treatment=="binary")
          else{
            bal[i,] <- abs(bal.tab(intervention~confounders,
                                   weights=gps_Slw$gps_w,
                                   method="weighting")$Balance$Corr.Adj)
          }#end else 
          
        }#end if NULL
      }#end if NULL
    }#end GPS_SL
    
    
    
    ####################################################################################################
    if(method=="BalSL"){
      if(type_treatment=="binary"){
        SL.library <- sl.lib_ps
        stop.method <- "es.mean"
        
        
        timet    <- tryCatch( system.time(ps_BalSL <- weightit(Tr ~ X.X1 + X.X2 + X.X3 + X.X4 + X.X5 + X.X6,
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
        
        timet    <- tryCatch( system.time(ps_BalSL <- weightit(Tr ~ X.X1 + X.X2 + X.X3 + X.X4 + X.X5 + X.X6,
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
          mse[i]  <- (coef[i]-true_haz)^2
          
          if(type_treatment=="binary"){
            bal[i,] <- abs(bal.tab(intervention~confounders,
                                   weights=ps_BalSL$weights,
                                   method="weighting",
                                   s.d.denom="pooled")$Balance$Diff.Adj)
          }#end if(type_treatment=="binary")
          else{
            bal[i,] <- abs(bal.tab(intervention~confounders,
                                   weights=ps_BalSL$weights,
                                   method="weighting")$Balance$Corr.Adj)
          }#end else 
          
        }#end if NULL
      }#end if NULL
    }#end Balance Super Learner 
    
    
    
    ####################################################################################################
    if(method=="GBM"){
      if(type_treatment=="binary"){
        stop.method <- "es.mean"
        
       
        timet    <- tryCatch( system.time(ps_gbm <- weightit(Tr ~ X.X1 + X.X2 + X.X3 + X.X4 + X.X5 + X.X6,
                                                               data = dta,
                                                               method = "gbm",
                                                               estimand = 'ATE',
                                                               stop.method = stop.method) )
                              , error=function(e) NULL)#endtryCatch
        
      }#end if treatment binary
      
      if(type_treatment=="continuous"){
        stop.method <- "p.mean"
        
        timet    <- tryCatch( system.time(ps_gbm <- weightit(Tr ~ X.X1 + X.X2 + X.X3 + X.X4 + X.X5 + X.X6,
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
          mse[i]  <- (coef[i]-true_haz)^2
          
          if(type_treatment=="binary"){
            bal[i,] <- abs(bal.tab(intervention~confounders,
                                   weights=ps_gbm$weights,
                                   method="weighting",
                                   s.d.denom="pooled")$Balance$Diff.Adj)
          }#end if(type_treatment=="binary")
          else{
            bal[i,] <- abs(bal.tab(intervention~confounders,
                                   weights=ps_gbm$weights,
                                   method="weighting")$Balance$Corr.Adj)
          }#end else 
          
        }#end if NULL
      }#end if NULL
    }#end GBM
    
    
    
    ####################################################################################################
    if(method=="SBW"){
      #SBW
      balsbw         <-  list()
      balsbw$bal_cov <- c("X.X1","X.X2","X.X3","X.X4","X.X5","X.X6")
      balsbw$bal_tol <- 0.01
      balsbw$bal_std <- TRUE
      
      timet <- tryCatch( system.time({capture.output( out0 <-  sbw(dta, ind = "Tr",  bal = balsbw,
                                                   sol = list(sol_nam = "gurobi"), par = list(par_est = "ate"))) 
                                      sbw <- out0$dat_weights$weights
                          })#end system.time 
      , error=function(e) NULL)#endtryCatch
      
      fit <- tryCatch(  coxph(Surv(time, status) ~ Tr , data=dta, weights=sbw)  
                 , error=function(e) NULL)#endtryCatch
      
      if(is.null(timet)!=TRUE){
          if(is.null(fit)!=TRUE){
            coef[i] <- fit$coefficients[1]
            mse[i]  <- (coef[i]-true_haz)^2
            #print(i)
            bal[i,] <- abs(bal.tab(intervention~confounders,
                                 weights=sbw,
                                 method="weighting",
                                 s.d.denom="pooled")$Balance$Diff.Adj)
        }#end if null fit
      }#end if null
      
    }#end method == "SBW"
    
    
    ####################################################################################################
    if(method == "PSM"){
      #compute ipw weight
      gps_Slw <- tryCatch( gps_weights(dta, type_treatment, sl.lib_gps, sl.lib_ps)
                            , error=function(e) NULL)#endtryCatch
      
      timet <- tryCatch(  system.time({ m1 <- Match(Tr = intervention, X = gps_Slw$gps_w, estimand = "ATE")
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
          mse[i]  <- (coef[i]-true_haz)^2
          bal[i,] <- abs(bal.tab(m1,intervention~confounders)$Balance$Diff.Adj)
          #bal[i,] <- abs(bal.tab(m2, disp.subclass = TRUE)$Balance$Diff.Adj[-1])
        }#end if null fit
        
      }#end if null

      
    }#end if method == PSM
    
    
    ####################################################################################################
    if(method == "EBAL"){
      
      #compute ipw weight
      timet <- tryCatch(  system.time({ebal   <- weightit(Tr ~ X.X1 +X.X2 +X.X3 +X.X4 +X.X5 +X.X6 , 
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
        mse[i]  <- (coef[i]-true_haz)^2
        bal[i,] <- abs(bal.tab(intervention~confounders,
                               weights=ebal$w,
                               method="weighting",
                               s.d.denom="pooled")$Balance$Diff.Adj)
        }#end if null
      }#end if null
    
    }#end if method == EBAL
    
    
    ####################################################################################################
    if(method=="OM"){
      timet        <- tryCatch( system.time(fit       <- summary(coxph(Surv(time, status) ~ Tr + X.X1 + X.X2 + X.X3 + X.X4 + X.X5 + X.X6, data=dta)))
                                , error=function(e) NULL)#endtryCatch

      if(is.null(timet)!=TRUE){
        coef[i]   <- fit$coefficients[1]
        mse[i]    <- (coef[i]-true_haz)^2
        bal[i,]   <- rep(0,6)
      }#end if NULL
    }#end OM
    
    
    ####################################################################################################
    if(method=="naive"){
      #naive
      timet      <- tryCatch( system.time(fit     <- summary(coxph(Surv(time, status) ~ Tr, data=dta)))
                              , error=function(e) NULL)#endtryCatch

      if(is.null(timet)!=TRUE){
        coef[i] <- fit$coefficients[1]
        mse[i]  <- (coef[i]-true_haz)^2
        
        if(type_treatment=="binary"){
          bal[i,] <- abs(bal.tab(intervention~confounders, s.d.denom = 'pooled')$Balance$Diff.Un)
        }#end if(type_treatment=="binary")
        else{
          bal[i,] <- abs(bal.tab(intervention~confounders)$Balance$Corr.Un)
        }#end else 
        

      }#end if NULL
    }#end method=="naive
    
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
  
}#end compute_haz_hat




compute_haz_hat_cova <- function(itera,method,rateCens,type_treatment,beta_pos,beta_miss,n,scenario,delta,num_cova){
  
  coef  <- mse <- time <- mcens <- timet <- NA
  bal   <- matrix(rep(NA,itera*num_cova),ncol=num_cova,nrow=itera)
  
  for(i in 1:itera){
    #print(i)
    if(i==1 | i%%25==0){print(paste("Iteration",i,"out of",itera))}  
    dta <- read.csv(paste("~/path",
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
          mse[i]  <- (coef[i]-true_haz)^2
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
                                  resampled_sample$X.X5,resampled_sample$X.X6)
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



#***********************************************************************
#
# This function provides standard errors as described in the Practical guidelines section
#
#***********************************************************************
compute_errors <- function(itera,method,rateCens,type_treatment,beta_pos,beta_miss,n,scenario,itera_boot,true_haz_exp,delta){
  
  sen <- ser <- seb <- coefn <- NA
  timeR <- timeN <- timeB <- NA
  cov95r <- cov95b <- cov95n <- NA
  mcens <- NA
  
  for(i in 1:itera){
    #print(i)
    if(i==1 | i%%25==0){print(paste("Iteration",i,"out of",itera))}  
    dta <- read.csv(paste("~/path",
                          rateCens,"_type_treatment_",type_treatment,"_lackover_",beta_pos,"_misspe_",beta_miss,"_scenario_",
                          scenario,"_n_",n,"_",i,".csv",sep=""))
    
    confounders   <- data.frame(dta$X.X1,dta$X.X2,dta$X.X3,dta$X.X4,dta$X.X5,dta$X.X6)
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
              cov95r[i] <- (tempS$conf.int[3]<=true_haz_exp)*(true_haz_exp<=tempS$conf.int[4])
            }#end if NULL
          })#end system timer
          
          ############
          #naive
          timen <- system.time({
            fitn          <- tryCatch( coxph(Surv(time, status) ~ Tr, data=dta, weights=row_w) , error=function(e) NULL)#endtryCatch
          
            if(is.null(fitn)!=TRUE){
              coefn[i]  <- fitn$coefficients[1]
              sen[i]    <- sqrt(fitn$var)
              tempS     <- summary(fitn)
              cov95n[i] <- (tempS$conf.int[3]<=true_haz_exp)*(true_haz_exp<=tempS$conf.int[4])
            }#end if NULL
          })#end system timer
          
      }#end if NULL
      
      ############
      #bootstrap
      timeb <- system.time({
        myBootstrap <- boot(dta, boot_row, R=itera_boot)
        seb[i]      <- sd(myBootstrap$t)
        ci          <- boot.ci(myBootstrap,conf=0.95,type="norm")
        cov95b[i]   <- (ci$normal[2]<=true_haz_exp)*(true_haz_exp<=ci$normal[3])
      })#end system timer
      
      timeR[i] <- timer[1] +  timer[2]
      timeN[i] <- timen[1] +  timen[2]
      timeB[i] <- timeb[1] +  timeb[2]
      
      mcens[i]  <- mean(dta$censoring,na.rm=T)
      Mcens   <- mean(mcens,na.rm=T)
      
      
    }#end method row
  }#end for (itera)
  
  if(is.null(coefn)!=T & is.null(sen)!=T & is.null(ser)!=T & is.null(seb)!=T  ){
    return(list(ses=sd(coefn,na.rm=T),
                sen=mean(sen,na.rm=T),
                seb=mean(seb,na.rm=T),
                ser=mean(ser,na.rm=T),
                timeR=mean(timeR,na.rm=T),
                timeN=mean(timeN,na.rm=T),
                timeB=mean(timeB,na.rm=T),
                cov95r=mean(cov95r,na.rm=T),
                cov95b=mean(cov95b,na.rm=T),
                cov95n=mean(cov95n,na.rm=T),
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








