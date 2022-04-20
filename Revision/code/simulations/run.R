########################################################################
########################################################################
#***********************************************************************
#
#   In this file, we provide functions for evaluating ROW across
#   levels of practical positivity violation, misspecification, and
#   censoring. In addition, we run simulations for coverage and standard 
#   errors
# 
#***********************************************************************
########################################################################
########################################################################


rm(list=ls())

set.seed(1)

library(ggpubr)
library(beepr)



###########Data generation variables
n             <- 1000
true_haz      <- 1.5
itera         <- 100
itera_boot    <- 100
depCensor     <- FALSE
depTreat      <- TRUE
seq_rateC     <- c(0.1,10)
#seq_miss      <- seq(0,1,length.out=2)
seq_miss      <- 1
seq_pos       <- c(1,2)
seq_scenario  <- c(1,2)
seq_type_trt  <- c("binary", "continuous")
#seq_type_trt  <- c("binary")
generate_data <- TRUE
delta         <- 1e-5


Bias <- MSE <- Bal <- Time <- Pcens <- Pscen <- Error <-  NULL
SES <- SEN <- SEB <- SEBW <- SER <-  TimeR <- TimeN <- TimeB <- TimeBW <- Co95B <- Co95BW <- Co95R <- Co95N <- NULL

sl.lib_gps <- c('SL.glm','SL.glm.interaction','SL.glm.BS','SL.glm.FP',
                'SL.glmnet',  'SL.rpart',"SL.bayesglm",
                "SL.loess",'SL.xgboost')
sl.lib_ps  <- c('SL.glm','SL.glm.interaction','SL.glm.FP',
                'SL.glmnet', 'SL.glm.BS', 'SL.rpart',"SL.bayesglm",'SL.xgboost')
ps          <- matrix(rep(NA,itera*n),nrow=itera,ncol=n)

Mbaln <- Mbalnc <- Mbalc <- Mbalw <- Mbalsl <- matrix(rep(NA,6*length(seq_rateC)),ncol=6,nrow=length(seq_rateC))

source("/Users/santam13/Documents/NYU/ROW/ROW-time-to-event/code/simulations/functions.R")









###################################################################################################################################################
###################################################################################################################################################
#
#
# Data Generation for simulation Hazard ratio and errors
#
#
###################################################################################################################################################
###################################################################################################################################################

if(isTRUE(generate_data)){
  
  
  ###################### **************** Positivity
  rateCens  <- 0.1  #Minimal Censoring
  beta_miss  <- 0    #No model misspecification
  scenario    <- 1 #Continuous treatment with minimal lack of overlap
  for(type_treatment in seq_type_trt){
    if(type_treatment=="binary"){
      for(beta_pos in seq_pos){
        print(paste(" ###################### **************** Lack of overlap Level ", beta_pos, " ****************  ###################### "))
        
        for(i in 1:itera){
          dta <- data_generation(n=n, true_haz = true_haz, 
                                 beta_tr=beta_pos, beta_miss=beta_miss,
                                 rateC=rateCens, 
                                 type_treatment=type_treatment, scenario=scenario)
          write.csv(dta, file=paste("/Users/santam13/Documents/NYU/ROW/ROW-time-to-event/code/simulations/data_tmp/",
                                    rateCens,"_type_treatment_",type_treatment,"_lackover_",beta_pos,"_misspe_",beta_miss,"_scenario_",
                                    scenario,"_n_",n,"_",i,".csv",sep=""), row.names=F)
        }#end for(i in 1:itera)
      }#end for(beta_pos in seq_pos){
    }#end if(type_treatment=="binary")
    else{
      beta_pos <- 1
      for(scenario in seq_scenario){
        print(paste(" ###################### **************** Lack of overlap Level ", scenario, " ****************  ###################### "))
        
        for(i in 1:itera){
          dta <- data_generation(n=n, true_haz = true_haz, 
                                 beta_tr=scenario, beta_miss=beta_miss,
                                 rateC=rateCens, 
                                 type_treatment=type_treatment, scenario=scenario)
          write.csv(dta, file=paste("/Users/santam13/Documents/NYU/ROW/ROW-time-to-event/code/simulations/data_tmp/",
                                    rateCens,"_type_treatment_",type_treatment,"_lackover_",beta_pos,"_misspe_",beta_miss,"_scenario_",
                                    scenario,"_n_",n,"_",i,".csv",sep=""), row.names=F)
        }#end for(i in 1:itera)
      }#end for(scenario in 1:3)
    }#end else if(type_treatment=="binary")
    
  }#end seq_type_trt
  
  
  ###################### **************** Misspecified 
  rateCens  <- 0.1 #Minimal Censoring
  beta_pos  <- 1.5 #Some lack of overlap for binary treatment
  scenario  <- 3 #Some lack of overlap for continuous treatment
  for(type_treatment in seq_type_trt){
    for(beta_miss in seq_miss){
      print(paste(" ###################### **************** Misspecification Level ", beta_miss, " ****************  ###################### "))
      
      for(i in 1:itera){
        dta <- data_generation(n=n, true_haz = true_haz, 
                               beta_tr=beta_pos, beta_miss=beta_miss,
                               rateC=rateCens, 
                               type_treatment=type_treatment, scenario=scenario)
        write.csv(dta, file=paste("/Users/santam13/Documents/NYU/ROW/ROW-time-to-event/code/simulations/data_tmp/",
                                  rateCens,"_type_treatment_",type_treatment,"_lackover_",beta_pos,"_misspe_",beta_miss,"_scenario_",
                                  scenario,"_n_",n,"_",i,".csv",sep=""), row.names=F)
      }#end for(i in 1:itera)
    }#end for(beta_miss in seq_miss){
  }#end seq_type_trt
  
  ###################### **************** Censoring 
  beta_miss   <- 0
  beta_pos    <- 1
  scenario    <- 1 #Continuous treatment with minimal lack of overlap
  for(type_treatment in seq_type_trt){
    for(rateCens in seq_rateC){
      print(paste(" ###################### **************** Rate Censoring ", rateCens, " ****************  ###################### "))
      
      for(i in 1:itera){
        dta <- data_generation(n=n, true_haz = true_haz, 
                               beta_tr=beta_pos, beta_miss=beta_miss,
                               rateC=rateCens, 
                               type_treatment=type_treatment, scenario=scenario)
        write.csv(dta, file=paste("/Users/santam13/Documents/NYU/ROW/ROW-time-to-event/code/simulations/data_tmp/",
                                  rateCens,"_type_treatment_",type_treatment,"_lackover_",beta_pos,"_misspe_",beta_miss,"_scenario_",
                                  scenario,"_n_",n,"_",i,".csv",sep=""), row.names=F)
      }#end for(i in 1:itera)
      
    }#end for(rateCens in seq_rateC)
  }#end seq_type_trt
  
}#end if(isTRUE(generate_data)){






#**************************************************************************************************************************************************
###################################################################################################################################################
###################################################################################################################################################
#
#
# Simulation for all methods
#
#
###################################################################################################################################################
###################################################################################################################################################
#**************************************************************************************************************************************************





###################################################################################################################################################
###################################################################################################################################################
#
#
# Simulation HAZARD RATIO
#
#
###################################################################################################################################################
###################################################################################################################################################

      seq_method_con <- c("IPW-SL","IPW-SL-stab","IPW-SL-trunc","BalSL","npCBPS","GBM","CBPS","IPW-BS","IPW-FP","OM","naive","ROW")
      seq_method_bin <- c("IPW-SL","IPW-SL-stab","IPW-SL-trunc","IPW-BS","IPW-FP","BalSL","SBW","GBM","PSM","CBPS","EBAL","OM","naive","ROW")
      # seq_type_trt <- c("binary", "continuous")
      # seq_type_trt <- c("binary")
      #
      # seq_method_bin <- c("IPW-BS")
      # seq_method_bin <- c("SBW","GBM", "ROW")
      # seq_method_con <- c("GBM", "ROW")
      # seq_type_trt <- c("continuous")
      # seq_type_trt <- c("binary")
      # seq_pos <- c(0.1,4)
      # sl.lib_gps <- c('SL.glm','SL.glm.interaction', 'SL.glm.BS', 'SL.rpart')
      # sl.lib_ps  <- c('SL.glm','SL.glm.interaction', 'SL.glm.BS', 'SL.rpart')
      
      
      ###################################################################################################################################################
      ###################################################################################################################################################
      # Positivity
      ###################################################################################################################################################
      ###################################################################################################################################################
      
      dfBIAS <- dfMSE <- dfBAL <- dfTIME <- dfERROR <- data.frame()
      
      name_bias <- c("Method", "Treatment", "Positivity", "Abs-Bias")
      name_mse  <- c("Method", "Treatment", "Positivity", "RMSE")
      name_bal  <- c("Method", "Treatment", "Positivity", "Balance")
      name_time <- c("Method", "Treatment", "Positivity", "Time")
      name_err  <- c("Method", "Treatment", "Positivity", "Error")
      
      rateCens  <- 0.1  #Minimal Censoring
      beta_miss  <- 0    #No model misspecification
      scenario    <- 1 #Continuous treatment with minimal lack of overlap
      for(type_treatment in seq_type_trt){
        gc()
        if(type_treatment=="binary"){
          seq_method <- seq_method_bin
        }#end if type_treatment
        else{
          seq_method <- seq_method_con
        }#end else
        
        for(method in seq_method){
        print(paste("########################### Method --->", method))
        if(type_treatment=="binary"){
          j <- 0
          for(beta_pos in seq_pos){
            j <- j+1
            print(paste(":::::::::::::::::::::::: Positivity",beta_pos, "- Treatment", type_treatment))
            res     <- compute_haz_hat(itera,method,rateCens,type_treatment,beta_pos,beta_miss,n,scenario,delta)
            Bias[j] <- res$bias
            MSE[j]  <- res$mse
            Bal[j]  <- mean(res$bal,na.rm=T)
            Time[j] <- res$time
            Pcens[j]<- res$pcens
            Error[j]<- mean(res$errors,na.rm=T)
            
          }#end for(beta_pos in seq_pos){
          
          
          ########################################################################
          #BIAS
          temp_bias           <- data.frame(rep(method,length(seq_pos)),
                                            rep(type_treatment,length(seq_pos)),
                                            seq_pos,Bias)
          colnames(temp_bias) <- name_bias
          dfBIAS              <- rbind(dfBIAS,temp_bias)
          
          ########################################################################
          #MSE
          temp_mse            <- data.frame(rep(method,length(seq_pos)),
                                            rep(type_treatment,length(seq_pos)),
                                            seq_pos,MSE)
          colnames(temp_mse)  <- name_mse
          dfMSE               <- rbind(dfMSE,temp_mse)
          
          ########################################################################
          #BAL
          temp_bal            <- data.frame(rep(method,length(seq_pos)),
                                            rep(type_treatment,length(seq_pos)),
                                            seq_pos,Bal)
          colnames(temp_bal)  <- name_bal
          dfBAL               <- rbind(dfBAL,temp_bal)
          
          ########################################################################
          #TIME
          temp_time           <- data.frame(rep(method,length(seq_pos)),
                                            rep(type_treatment,length(seq_pos)),
                                            seq_pos,Time)
          colnames(temp_time) <- name_time
          dfTIME              <- rbind(dfTIME,temp_time)
          
          ########################################################################
          #ERROR
          temp_err            <- data.frame(rep(method,length(seq_pos)),
                                            rep(type_treatment,length(seq_pos)),
                                            seq_pos,Error)
          colnames(temp_err)  <- name_err
          dfERROR             <- rbind(dfERROR,temp_err)
          
        }#end if(type_treatment=="binary")
        else{
          j <- 0
          beta_pos <- 1
          for(scenario in seq_scenario){
            j <- j+1
            print(paste(":::::::::::::::::::::::: Positivity",scenario, "- Treatment", type_treatment))
            res     <- compute_haz_hat(itera,method,rateCens,type_treatment,beta_pos,beta_miss,n,scenario,delta)
            Bias[j] <- res$bias
            MSE[j]  <- res$mse
            Bal[j]  <- mean(res$bal,na.rm=T)
            Time[j] <- res$time
            Pcens[j]<- res$pcens
            Pscen[j]<- scenario
            Error[j]<- mean(res$errors,na.rm=T)
            
          }#end for(scenario in 1:3)
          
          ########################################################################
          #BIAS
          temp_bias           <- data.frame(rep(method,length(seq_scenario)),
                                            rep(type_treatment,length(seq_scenario)),
                                            Pscen,Bias)
          colnames(temp_bias) <- name_bias
          dfBIAS              <- rbind(dfBIAS,temp_bias)
          
          ########################################################################
          #MSE
          temp_mse            <- data.frame(rep(method,length(seq_scenario)),
                                            rep(type_treatment,length(seq_scenario)),
                                            Pscen,MSE)
          colnames(temp_mse)  <- name_mse
          dfMSE               <- rbind(dfMSE,temp_mse)
          
          ########################################################################
          #BAL
          temp_bal            <- data.frame(rep(method,length(seq_scenario)),
                                            rep(type_treatment,length(seq_scenario)),
                                            Pscen,Bal)
          colnames(temp_bal)  <- name_bal
          dfBAL               <- rbind(dfBAL,temp_bal)
          
          ########################################################################
          #TIME
          temp_time           <- data.frame(rep(method,length(seq_scenario)),
                                            rep(type_treatment,length(seq_scenario)),
                                            Pscen,Time)
          colnames(temp_time) <- name_time
          dfTIME              <- rbind(dfTIME,temp_time)
          
          ########################################################################
          #ERROR
          temp_err            <- data.frame(rep(method,length(seq_pos)),
                                            rep(type_treatment,length(seq_pos)),
                                            seq_pos,Error)
          colnames(temp_err)  <- name_err
          dfERROR             <- rbind(dfERROR,temp_err)
          
        }#end else if(type_treatment=="binary")
          save.image(paste("/Users/santam13/Documents/NYU/ROW/ROW-time-to-event/code/simulations/data_results/simu_n",n,"positivity_nov012021.Rdata"))
          
        }#end for(method in seq_method)
      }#end seq_type_trt
      
      
      save.image(paste("/Users/santam13/Documents/NYU/ROW/ROW-time-to-event/code/simulations/data_results/simu_n",n,"positivity.Rdata"))
      
      # beep(10)
      
      
      
      ###################################################################################################################################################
      ###################################################################################################################################################
      # Misspecification
      ###################################################################################################################################################
      ###################################################################################################################################################
      
      
      dfBIAS <- dfMSE <- dfBAL <- dfTIME <- dfERROR <- data.frame()
      
      name_bias <- c("Method", "Treatment", "Misspecification", "Abs-Bias")
      name_mse  <- c("Method", "Treatment", "Misspecification", "RMSE")
      name_bal  <- c("Method", "Treatment", "Misspecification", "Balance")
      name_time <- c("Method", "Treatment", "Misspecification", "Time")
      name_err  <- c("Method", "Treatment", "Misspecification", "Error")
      
      rateCens  <- 0.1 #Minimal Censoring
      beta_pos  <- 1.5 #Some lack of overlap for binary treatment
      scenario  <- 3 #Some lack of overlap for continuous treatment
      for(type_treatment in seq_type_trt){
        
        if(type_treatment=="binary"){
          seq_method <- seq_method_bin
        }#end if type_treatment
        else{
          seq_method <- seq_method_con
        }#end else
        
        for(method in seq_method){
          print(paste("########################### Method --->", method))
            j <- 0
            for(beta_miss in seq_miss){
              j <- j+1
              print(paste(":::::::::::::::::::::::: Misspecification",beta_miss, "- Treatment", type_treatment))
              res     <- compute_haz_hat(itera,method,rateCens,type_treatment,beta_pos,beta_miss,n,scenario,delta)
              Bias[j] <- res$bias
              MSE[j]  <- res$mse
              Bal[j]  <- mean(res$bal,na.rm=T)
              Time[j] <- res$time
              Pcens[j]<- res$pcens
              Error[j]<- mean(res$errors,na.rm=T)
      
            }#end for(beta_pos in seq_pos){
            
            
            ########################################################################
            #BIAS
            temp_bias           <- data.frame(rep(method,length(seq_miss)),
                                              rep(type_treatment,length(seq_miss)),
                                              seq_miss,Bias)
            colnames(temp_bias) <- name_bias
            dfBIAS              <- rbind(dfBIAS,temp_bias)
            
            ########################################################################
            #MSE
            temp_mse            <- data.frame(rep(method,length(seq_miss)),
                                              rep(type_treatment,length(seq_miss)),
                                              seq_miss,MSE)
            colnames(temp_mse)  <- name_mse
            dfMSE               <- rbind(dfMSE,temp_mse)
            
            ########################################################################
            #BAL
            temp_bal            <- data.frame(rep(method,length(seq_miss)),
                                              rep(type_treatment,length(seq_miss)),
                                              seq_miss,Bal)
            colnames(temp_bal)  <- name_bal
            dfBAL               <- rbind(dfBAL,temp_bal)
            
            ########################################################################
            #TIME
            temp_time           <- data.frame(rep(method,length(seq_miss)),
                                              rep(type_treatment,length(seq_miss)),
                                              seq_miss,Time)
            colnames(temp_time) <- name_time
            dfTIME              <- rbind(dfTIME,temp_time)
            
            ########################################################################
            #ERROR
            temp_err            <- data.frame(rep(method,length(seq_pos)),
                                              rep(type_treatment,length(seq_pos)),
                                              seq_pos,Error)
            colnames(temp_err)  <- name_err
            dfERROR             <- rbind(dfERROR,temp_err)
            
          
        }#end for(method in seq_method)
      }#end seq_type_trt
      
      
      
      save.image(paste("/Users/santam13/Documents/NYU/ROW/ROW-time-to-event/code/simulations/data_results/simu_n",n,"misspecified.Rdata"))
      
      # beep(9)
      
      
      
      
      ###################################################################################################################################################
      ###################################################################################################################################################
      # Censoring
      ###################################################################################################################################################
      ###################################################################################################################################################
      
      
      dfBIAS <- dfMSE <- dfBAL <- dfTIME <- data.frame()
      
      name_bias <- c("Method", "Treatment", "Censoring", "Bias")
      name_mse  <- c("Method", "Treatment", "Censoring", "MSE")
      name_bal  <- c("Method", "Treatment", "Censoring", "Balance")
      name_time <- c("Method", "Treatment", "Censoring", "Time")
      
      beta_miss   <- 0
      beta_pos    <- 1
      scenario    <- 1 #Continuous treatment with minimal lack of overlap
      for(type_treatment in seq_type_trt){
        
        if(type_treatment=="binary"){
          seq_method <- seq_method_bin
        }#end if type_treatment
        else{
          seq_method <- seq_method_con
        }#end else
        
        for(method in seq_method){
          print(paste("########################### Method --->", method))
          j <- 0
          for(rateCens in seq_rateC){
            j <- j+1
            print(paste(":::::::::::::::::::::::: Censoring",rateCens, "- Treatment", type_treatment))
            res     <- compute_haz_hat(itera,method,rateCens,type_treatment,beta_pos,beta_miss,n,scenario,delta)
            Bias[j] <- res$bias
            MSE[j]  <- res$mse
            Bal[j]  <- mean(res$bal,na.rm=T)
            Time[j] <- res$time
            Pcens[j]<- res$pcens
            
          }#end for(beta_pos in seq_pos){
          
          
          ########################################################################
          #BIAS
          temp_bias           <- data.frame(rep(method,length(seq_rateC)),
                                            rep(type_treatment,length(seq_rateC)),
                                            Pcens,Bias)
          colnames(temp_bias) <- name_bias
          dfBIAS              <- rbind(dfBIAS,temp_bias)
          
          ########################################################################
          #MSE
          temp_mse            <- data.frame(rep(method,length(seq_rateC)),
                                            rep(type_treatment,length(seq_rateC)),
                                            Pcens,MSE)
          colnames(temp_mse)  <- name_mse
          dfMSE               <- rbind(dfMSE,temp_mse)
          
          ########################################################################
          #BAL
          temp_bal            <- data.frame(rep(method,length(seq_rateC)),
                                            rep(type_treatment,length(seq_rateC)),
                                            Pcens,Bal)
          colnames(temp_bal)  <- name_bal
          dfBAL               <- rbind(dfBAL,temp_bal)
          
          ########################################################################
          #TIME
          temp_time           <- data.frame(rep(method,length(seq_rateC)),
                                            rep(type_treatment,length(seq_rateC)),
                                            Pcens,Time)
          colnames(temp_time) <- name_time
          dfTIME              <- rbind(dfTIME,temp_time)
          
          
        }#end for(method in seq_method)
      }#end seq_type_trt
      
      save.image(paste("/Users/santam13/Documents/NYU/ROW/ROW-time-to-event/code/simulations/data_results/simu_n",n,"censoring.Rdata"))



      



#**************************************************************************************************************************************************
###################################################################################################################################################
###################################################################################################################################################
#
#
# Simulation Only for ROW
#
#
###################################################################################################################################################
###################################################################################################################################################
#**************************************************************************************************************************************************
      





###################################################################################################################################################
###################################################################################################################################################
#
#
# Simulation COVERAGE - SE
#
#
###################################################################################################################################################
###################################################################################################################################################


      seq_method_bin <- seq_method_con <- "ROW"
      
      
      ###################################################################################################################################################
      ###################################################################################################################################################
      # Positivity
      ###################################################################################################################################################
      ###################################################################################################################################################
      
      SES <- SEN <- SEB <- SEBW <- SER <-  TimeR <- TimeN <- TimeB <- TimeBW <- Co95B <- Co95BW <- Co95R <- Co95N <- Co95B_n <- Co95BW_n <- NULL
      dfSE <- dfCO <- dfTIME <- data.frame()
      
      name_se   <- c("Type","Treatment", "Positivity", "SE")
      name_co   <- c("Type","Treatment", "Positivity", "Coverage")
      name_ti   <- c("Type","Treatment", "Positivity", "Time")
      
      delta <- 1e-05
      rateCens  <- 0.1  #Minimal Censoring
      beta_miss  <- 0    #No model misspecification
      scenario <- 1 #Minimal lack of overlap for continuous treatment
      true_haz_exp <- exp(true_haz)
      for(type_treatment in seq_type_trt){
        
        if(type_treatment=="binary"){
          seq_method <- seq_method_bin
        }#end if type_treatment
        else{
          seq_method <- seq_method_con
        }#end else
        
        for(method in seq_method){
          print(paste("########################### Method --->", method))
          if(type_treatment=="binary"){
            j <- 0
            for(beta_pos in seq_pos){
              j <- j+1
              print(paste(":::::::::::::::::::::::: Positivity",beta_pos, "- Treatment", type_treatment))
              res     <- compute_errors(itera,method,rateCens,type_treatment,beta_pos,beta_miss,n,scenario,itera_boot,true_haz,delta)
              SES[j]  <- res$ses
              SEB[j]  <- res$seb
              SEBW[j]  <- res$sebw
              SEN[j]  <- res$sen
              SER[j]  <- res$ser
              TimeR[j]<- res$timeR
              TimeN[j]<- res$timeN
              TimeB[j]<- res$timeB
              TimeBW[j]<- res$timeBW
              Co95B[j]<- res$cov95b
              Co95B_n[j]<- res$cov95b_n
              Co95BW[j]<- res$cov95bw
              Co95BW_n[j]<- res$cov95bw_n
              Co95R[j]<- res$cov95r
              Co95N[j]<- res$cov95n
              
              
            }#end for(beta_pos in seq_pos){
            
            
            ########################################################################
            #SE
            temp_se           <- data.frame( c(rep('Empirical',length(seq_pos)),
                                               rep('Bootstrap',length(seq_pos)),
                                               rep('W-Bootstrap',length(seq_pos)),
                                               rep('Robust',length(seq_pos)),
                                               rep('Naive',length(seq_pos) ) 
                                               ),
                                              rep(type_treatment,length(seq_pos)),
                                              rep(seq_pos,5),
                                              c(SES,SEB,SEBW,SER,SEN))
            colnames(temp_se) <- name_se
            dfSE            <- rbind(dfSE,temp_se)
            
            ########################################################################
            #Coverage
            temp_co           <- data.frame( c(rep('Bootstrap-Perc',length(seq_pos)),
                                               rep('Bootstrap-Norm',length(seq_pos)),
                                               rep('W-Bootstrap-Perc',length(seq_pos)),
                                               rep('W-Bootstrap-Norm',length(seq_pos)),
                                               rep('Robust',length(seq_pos)),
                                               rep('Naive',length(seq_pos) ) 
                                              ),
                                              rep(type_treatment,length(seq_pos)),
                                              rep(seq_pos,6),
                                              c(Co95B,Co95B_n,Co95BW,Co95BW_n,Co95R,Co95N))
            colnames(temp_co) <- name_co
            dfCO            <- rbind(dfCO,temp_co)
            
            ########################################################################
            #TIME
            temp_ti           <- data.frame( c(rep('Bootstrap',length(seq_pos)),
                                               rep('W-Bootstrap',length(seq_pos)),
                                               rep('Robust',length(seq_pos)),
                                               rep('Naive',length(seq_pos) ) 
                                            ),
                                            rep(type_treatment,length(seq_pos)),
                                            rep(seq_pos,4),
                                            c(TimeB,TimeBW,TimeR,TimeN))
            colnames(temp_ti) <- name_ti
            dfTIME            <- rbind(dfTIME,temp_ti)
            
          }#end if(type_treatment=="binary")
          else{
            j <- 0
            beta_pos <- 1
            for(scenario in seq_scenario){
              j <- j+1
              print(paste(":::::::::::::::::::::::: Positivity",scenario, "- Treatment", type_treatment))
              res     <- compute_errors(itera,method,rateCens,type_treatment,beta_pos,beta_miss,n,scenario,itera_boot,true_haz,delta)
              SES[j]  <- res$ses
              SEB[j]  <- res$seb
              SEBW[j]  <- res$sebw
              SEN[j]  <- res$sen
              SER[j]  <- res$ser
              TimeR[j]<- res$timeR
              TimeN[j]<- res$timeN
              TimeB[j]<- res$timeB
              TimeBW[j]<- res$timeBW
              Co95B[j]<- res$cov95b
              Co95B_n[j]<- res$cov95b_n
              Co95BW[j]<- res$cov95bw
              Co95BW_n[j]<- res$cov95bw_n
              Co95R[j]<- res$cov95r
              Co95N[j]<- res$cov95n
              
            }#end for(scenario in 1:3)
            
            
            ########################################################################
            #SE
            temp_se           <- data.frame( c(rep('Empirical',length(seq_pos)),
                                               rep('Bootstrap',length(seq_pos)),
                                               rep('W-Bootstrap',length(seq_pos)),
                                               rep('Robust',length(seq_pos)),
                                               rep('Naive',length(seq_pos) ) 
            ),
            rep(type_treatment,length(seq_pos)),
            rep(seq_pos,5),
            c(SES,SEB,SEBW,SER,SEN))
            colnames(temp_se) <- name_se
            dfSE            <- rbind(dfSE,temp_se)
            
            ########################################################################
            #Coverage
            temp_co           <- data.frame( c(rep('Bootstrap-Perc',length(seq_pos)),
                                               rep('Bootstrap-Norm',length(seq_pos)),
                                               rep('W-Bootstrap-Perc',length(seq_pos)),
                                               rep('W-Bootstrap-Norm',length(seq_pos)),
                                               rep('Robust',length(seq_pos)),
                                               rep('Naive',length(seq_pos) ) 
            ),
            rep(type_treatment,length(seq_pos)),
            rep(seq_pos,6),
            c(Co95B,Co95B_n,Co95BW,Co95BW_n,Co95R,Co95N))
            colnames(temp_co) <- name_co
            dfCO            <- rbind(dfCO,temp_co)
            
            ########################################################################
            #TIME
            temp_ti           <- data.frame( c(rep('Bootstrap',length(seq_pos)),
                                               rep('W-Bootstrap',length(seq_pos)),
                                               rep('Robust',length(seq_pos)),
                                               rep('Naive',length(seq_pos) ) 
            ),
            rep(type_treatment,length(seq_pos)),
            rep(seq_pos,4),
            c(TimeB,TimeBW,TimeR,TimeN))
            colnames(temp_ti) <- name_ti
            dfTIME            <- rbind(dfTIME,temp_ti)
          }#end else if(type_treatment=="binary")
          
        }#end for(method in seq_method)
      }#end seq_type_trt
      
      save.image(paste("/Users/santam13/Documents/NYU/ROW/ROW-time-to-event/code/simulations/data_results/simu_errors_n",n,"positivity.Rdata"))
      
      
      
      ###################################################################################################################################################
      ###################################################################################################################################################
      # Misspecification
      ###################################################################################################################################################
      ###################################################################################################################################################
      
      SES <- SEN <- SEB <- SEBW <- SER <-  TimeR <- TimeN <- TimeB <- TimeBW <- Co95B <- Co95BW <- Co95R <- Co95N <- Co95B_n <- Co95BW_n <- NULL
      dfSE <- dfCO <- dfTIME <- data.frame()
      
      name_se   <- c("Type","Treatment", "Misspecification", "SE")
      name_co   <- c("Type","Treatment", "Misspecification", "Coverage")
      name_ti   <- c("Type","Treatment", "Misspecification", "Time")
      
      
      rateCens  <- 0.1 #Minimal Censoring
      beta_pos  <-  1.5 #Moderate lack of overlap for binary treatment
      scenario    <- 3 #Moderate lack of overlap for continuous treatment
      true_haz_exp <- exp(true_haz)
      for(type_treatment in seq_type_trt){
        
        if(type_treatment=="binary"){
          seq_method <- seq_method_bin
        }#end if type_treatment
        else{
          seq_method <- seq_method_con
        }#end else
        
        for(method in seq_method){
          print(paste("########################### Method --->", method))
            j <- 0
            for(beta_miss in seq_miss){
              j <- j+1
              print(paste(":::::::::::::::::::::::: Misspecification",beta_miss, "- Treatment", type_treatment))
              res     <- compute_errors(itera,method,rateCens,type_treatment,beta_pos,beta_miss,n,scenario,itera_boot,true_haz,delta)
              SES[j]  <- res$ses
              SEB[j]  <- res$seb
              SEBW[j]  <- res$sebw
              SEN[j]  <- res$sen
              SER[j]  <- res$ser
              TimeR[j]<- res$timeR
              TimeN[j]<- res$timeN
              TimeB[j]<- res$timeB
              TimeBW[j]<- res$timeBW
              Co95B[j]<- res$cov95b
              Co95B_n[j]<- res$cov95b_n
              Co95BW[j]<- res$cov95bw
              Co95BW_n[j]<- res$cov95bw_n
              Co95R[j]<- res$cov95r
              Co95N[j]<- res$cov95n
              
              
            }#end for(beta_pos in seq_miss){
            
            
            ########################################################################
            #SE
            temp_se           <- data.frame( c(rep('Empirical',length(seq_miss)),
                                               rep('Bootstrap',length(seq_miss)),
                                               rep('W-Bootstrap',length(seq_miss)),
                                               rep('Robust',length(seq_miss)),
                                               rep('Naive',length(seq_miss) ) 
            ),
            rep(type_treatment,length(seq_miss)),
            rep(seq_miss,5),
            c(SES,SEB,SEBW,SER,SEN))
            colnames(temp_se) <- name_se
            dfSE            <- rbind(dfSE,temp_se)
            
            ########################################################################
            #Coverage
            temp_co           <- data.frame( c(rep('Bootstrap-Perc',length(seq_miss)),
                                               rep('Bootstrap-Norm',length(seq_miss)),
                                               rep('W-Bootstrap-Perc',length(seq_miss)),
                                               rep('W-Bootstrap-Norm',length(seq_miss)),
                                               rep('Robust',length(seq_miss)),
                                               rep('Naive',length(seq_miss) ) 
            ),
            rep(type_treatment,length(seq_miss)),
            rep(seq_miss,6),
            c(Co95B,Co95B_n,Co95BW,Co95BW_n,Co95R,Co95N))
            colnames(temp_co) <- name_co
            dfCO            <- rbind(dfCO,temp_co)
            
            ########################################################################
            #TIME
            temp_ti           <- data.frame( c(rep('Bootstrap',length(seq_miss)),
                                               rep('W-Bootstrap',length(seq_miss)),
                                               rep('Robust',length(seq_miss)),
                                               rep('Naive',length(seq_miss) ) 
            ),
            rep(type_treatment,length(seq_miss)),
            rep(seq_miss,4),
            c(TimeB,TimeBW,TimeR,TimeN))
            colnames(temp_ti) <- name_ti
            dfTIME            <- rbind(dfTIME,temp_ti)
        
          
        }#end for(method in seq_method)
      }#end seq_type_trt
      
      save.image(paste("/Users/santam13/Documents/NYU/ROW/ROW-time-to-event/code/simulations/data_results/simu_errors_n",n,"misspecification.Rdata"))
      
      

      ###################################################################################################################################################
      ###################################################################################################################################################
      # Censoring
      ###################################################################################################################################################
      ###################################################################################################################################################
      
      SES <- SEN <- SEB <- SEBW <- SER <-  TimeR <- TimeN <- TimeB <- TimeBW <- Co95B <- Co95BW <- Co95R <- Co95N <- NULL
      dfSE <- dfCO <- dfTIME <- data.frame()
      
      name_se   <- c("Type","Treatment", "Censoring", "SE")
      name_co   <- c("Type","Treatment", "Censoring", "Coverage")
      name_ti   <- c("Type","Treatment", "Censoring", "Time")
      
      
      beta_miss  <- 0    #No model misspecification
      beta_pos  <- 1 #Moderate lack of overlap for binary treatment
      scenario    <- 1 #Minimal lack of overlap for continuous treatment
      true_haz_exp <- exp(true_haz)
      for(type_treatment in seq_type_trt){
        
        if(type_treatment=="binary"){
          seq_method <- seq_method_bin
        }#end if type_treatment
        else{
          seq_method <- seq_method_con
        }#end else
        
        for(method in seq_method){
          print(paste("########################### Method --->", method))
          j <- 0
          for(rateCens in seq_rateC){
            j <- j+1
            print(paste(":::::::::::::::::::::::: Censoring",rateCens, "- Treatment", type_treatment))
            res     <- compute_errors(itera,method,rateCens,type_treatment,beta_pos,beta_miss,n,scenario,itera_boot,true_haz,delta)
            SES[j]  <- res$ses
            SEB[j]  <- res$seb
            SEBW[j]  <- res$sebw
            SEN[j]  <- res$sen
            SER[j]  <- res$ser
            TimeR[j]<- res$timeR
            TimeN[j]<- res$timeN
            TimeB[j]<- res$timeB
            TimeBW[j]<- res$timeBW
            Co95B[j]<- res$cov95b
            Co95BW[j]<- res$cov95bw
            Co95R[j]<- res$cov95r
            Co95N[j]<- res$cov95n
            
          }#end for(beta_pos in seq_rateC){
          
          
          ########################################################################
          #SE
          temp_se           <- data.frame( c(rep('Empirical',length(seq_rateC)),
                                             rep('Bootstrap',length(seq_rateC)),
                                             rep('W-Bootstrap',length(seq_rateC)),
                                             rep('Robust',length(seq_rateC)),
                                             rep('Naive',length(seq_rateC) ) 
          ),
          rep(type_treatment,length(seq_rateC)),
          rep(seq_rateC,5),
          c(SES,SEB,SEBW,SER,SEN))
          colnames(temp_se) <- name_se
          dfSE            <- rbind(dfSE,temp_se)
          
          ########################################################################
          #Coverage
          temp_co           <- data.frame( c(rep('Bootstrap',length(seq_rateC)),
                                             rep('W-Bootstrap',length(seq_rateC)),
                                             rep('Robust',length(seq_rateC)),
                                             rep('Naive',length(seq_rateC) ) 
          ),
          rep(type_treatment,length(seq_rateC)),
          rep(seq_rateC,4),
          c(Co95B,Co95BW,Co95R,Co95N))
          colnames(temp_co) <- name_co
          dfCO            <- rbind(dfCO,temp_co)
          
          ########################################################################
          #TIME
          temp_ti           <- data.frame( c(rep('Bootstrap',length(seq_rateC)),
                                             rep('W-Bootstrap',length(seq_rateC)),
                                             rep('Robust',length(seq_rateC)),
                                             rep('Naive',length(seq_rateC) ) 
          ),
          rep(type_treatment,length(seq_rateC)),
          rep(seq_rateC,4),
          c(TimeB,TimeBW,TimeR,TimeN))
          colnames(temp_ti) <- name_ti
          dfTIME            <- rbind(dfTIME,temp_ti)
          
          
          
        }#end for(method in seq_method)
      }#end seq_type_trt
      
      save.image(paste("/Users/santam13/Documents/NYU/ROW/ROW-time-to-event/code/simulations/data_results/simu_errors_n",n,"censoring.Rdata"))
      
      
      
  ###################################################################################################################################################
  ###################################################################################################################################################
  #
  #
  # Simulation OVER delta
  #
  #
  ###################################################################################################################################################
  ###################################################################################################################################################
  
      n         <- 1000
      seq_method<- "ROW"
      #seq_delta <- c(1e-08, 1e-05, 1e-04, 1e-03, 1e-02, 1e-01, 0.2, 0.5) 
      seq_delta <- c(1e-08, 1e-05, 0.025, 0.05, 0.075, 0.1, 0.2, 0.5) 
      
      ###################################################################################################################################################
      ###################################################################################################################################################
      # Simulations over delta
      ###################################################################################################################################################
      ###################################################################################################################################################
      
      dfBIAS <- dfMSE <- dfBAL <- dfTIME <- dfERROR <- data.frame()
      
      name_bias <- c("Method", "Treatment", "Delta", "Bias")
      name_mse  <- c("Method", "Treatment", "Delta", "MSE")
      name_bal  <- c("Method", "Treatment", "Delta", "Balance")
      name_time <- c("Method", "Treatment", "Delta", "Time")
      name_err  <- c("Method", "Treatment", "Delta", "Error")
      
      num_cova <- 4
      #USE Data generation over n
      
      rateCens  <- 0.1  #Minimal Censoring
      beta_miss  <- 0    #No model misspecification
      beta_pos <- 1
      scenario    <- 1 #Continuous treatment with minimal lack of overlap
      for(type_treatment in seq_type_trt){
        
        if(type_treatment=="binary"){
          seq_method <- seq_method_bin
        }#end if type_treatment
        else{
          seq_method <- seq_method_con
        }#end else
        
        for(method in seq_method){
          print(paste("########################### Method --->", method))
          j <- 0
          for(delta in seq_delta){
            j <- j+1
            print(paste(":::::::::::::::::::::::: Delta", delta, "- Treatment", type_treatment))
            #res     <- compute_haz_hat_cova(itera,method,rateCens,type_treatment,beta_pos,beta_miss,n,scenario,delta,num_cova)
            res <- compute_haz_hat(itera,method,rateCens,type_treatment,beta_pos,beta_miss,n,scenario,delta)
            Bias[j] <- res$bias
            MSE[j]  <- res$mse
            Bal[j]  <- mean(res$bal,na.rm=T)
            Time[j] <- res$time
            Pcens[j]<- res$pcens
            Error[j]<- mean(res$errors,na.rm=T)
            
            
          }#end for(beta_pos in seq_pos){
          
          
          ########################################################################
          #BIAS
          temp_bias           <- data.frame(rep(method,length(seq_delta)),
                                            rep(type_treatment,length(seq_delta)),
                                            seq_delta,Bias)
          colnames(temp_bias) <- name_bias
          dfBIAS              <- rbind(dfBIAS,temp_bias)
          
          ########################################################################
          #MSE
          temp_mse            <- data.frame(rep(method,length(seq_delta)),
                                            rep(type_treatment,length(seq_delta)),
                                            seq_delta,MSE)
          colnames(temp_mse)  <- name_mse
          dfMSE               <- rbind(dfMSE,temp_mse)
          
          ########################################################################
          #BAL
          temp_bal            <- data.frame(rep(method,length(seq_delta)),
                                            rep(type_treatment,length(seq_delta)),
                                            seq_delta,Bal)
          colnames(temp_bal)  <- name_bal
          dfBAL               <- rbind(dfBAL,temp_bal)
          
          ########################################################################
          #TIME
          temp_time           <- data.frame(rep(method,length(seq_delta)),
                                            rep(type_treatment,length(seq_delta)),
                                            seq_delta,Time)
          colnames(temp_time) <- name_time
          dfTIME              <- rbind(dfTIME,temp_time)
          
          
          ########################################################################
          #ERROR
          temp_err            <- data.frame(rep(method,length(seq_pos)),
                                            rep(type_treatment,length(seq_pos)),
                                            seq_pos,Error)
          colnames(temp_err)  <- name_err
          dfERROR             <- rbind(dfERROR,temp_err)
          
          
        }#end for(method in seq_method)
      }#end seq_type_trt
      
      save.image(paste("/Users/santam13/Documents/NYU/ROW/ROW-time-to-event/code/simulations/data_results/simu_n",n,"delta.Rdata"))
          
      
      
      
      
      
      

###################################################################################################################################################
###################################################################################################################################################
#
#
# Simulation OVER N and OVER COVARIATES
#
#
###################################################################################################################################################
###################################################################################################################################################

      set.seed(1)
      seq_nn <- c(50,75,100,250,500,1000,5000,10000)
      seq_method_bin <- seq_method_con <- "ROW"
      
      
      num_cova <- 4
      
      ###################################################################################################################################################
      ###################################################################################################################################################
      # Data generation over n
      ###################################################################################################################################################
      ###################################################################################################################################################
      
      ###################### **************** Over n
      beta_miss   <- 0
      beta_pos    <- 1
      scenario    <- 1 #Continuous treatment with minimal lack of overlap
      for(type_treatment in seq_type_trt){
        for(n in seq_nn){
          print(paste(" ###################### **************** N ", n, " ****************  ###################### "))
          
          for(i in 1:itera){
            dta <- data_generation_cova(n=n, lambda=lambda, rho=rho, delta = true_haz, 
                                        beta_tr=beta_pos, beta_miss=beta_miss,
                                        beta_ou=beta_ou, rateC=rateCens, 
                                        type_treatment=type_treatment, scenario=scenario,
                                        num_cova=num_cova)
            write.csv(dta, file=paste("~/path",
                                      rateCens,"_type_treatment_",type_treatment,"_lackover_",beta_pos,"_misspe_",beta_miss,"_scenario_",
                                      scenario,"_n_",n,"_","num_cova",num_cova,"_",i,".csv",sep=""), row.names=F)
          }#end for(i in 1:itera)
          
        }#end for(rateCens in seq_rateC)
      }#end seq_type_trt
      
      
      
      ###################################################################################################################################################
      ###################################################################################################################################################
      # Simulations over n
      ###################################################################################################################################################
      ###################################################################################################################################################
      
      dfBIAS <- dfMSE <- dfBAL <- dfTIME <- data.frame()
      
      name_bias <- c("Method", "Treatment", "Sample size", "Bias")
      name_mse  <- c("Method", "Treatment", "Sample size", "MSE")
      name_bal  <- c("Method", "Treatment", "Sample size", "Balance")
      name_time <- c("Method", "Treatment", "Sample size", "Time")
      
      beta_miss  <- 0    #No model misspecification
      beta_pos  <- 1 #Moderate lack of overlap for binary treatment
      scenario    <- 1 #Minimal lack of overlap for continuous treatment
      for(type_treatment in seq_type_trt){
        
        if(type_treatment=="binary"){
          seq_method <- seq_method_bin
        }#end if type_treatment
        else{
          seq_method <- seq_method_con
        }#end else
        
        for(method in seq_method){
          print(paste("########################### Method --->", method))
          j <- 0
          for(n in seq_nn){
            j <- j+1
            print(paste(":::::::::::::::::::::::: N", n, "- Treatment", type_treatment))
            res     <- compute_haz_hat_cova(itera,method,rateCens,type_treatment,beta_pos,beta_miss,n,scenario,delta,num_cova)
            Bias[j] <- res$bias
            MSE[j]  <- res$mse
            Bal[j]  <- mean(res$bal,na.rm=T)
            Time[j] <- res$time
            Pcens[j]<- res$pcens
            
          }#end for(beta_pos in seq_pos){
          
          
          ########################################################################
          #BIAS
          temp_bias           <- data.frame(rep(method,length(seq_nn)),
                                            rep(type_treatment,length(seq_nn)),
                                            seq_nn,Bias)
          colnames(temp_bias) <- name_bias
          dfBIAS              <- rbind(dfBIAS,temp_bias)
          
          ########################################################################
          #MSE
          temp_mse            <- data.frame(rep(method,length(seq_nn)),
                                            rep(type_treatment,length(seq_nn)),
                                            seq_nn,MSE)
          colnames(temp_mse)  <- name_mse
          dfMSE               <- rbind(dfMSE,temp_mse)
          
          ########################################################################
          #BAL
          temp_bal            <- data.frame(rep(method,length(seq_nn)),
                                            rep(type_treatment,length(seq_nn)),
                                            seq_nn,Bal)
          colnames(temp_bal)  <- name_bal
          dfBAL               <- rbind(dfBAL,temp_bal)
          
          ########################################################################
          #TIME
          temp_time           <- data.frame(rep(method,length(seq_nn)),
                                            rep(type_treatment,length(seq_nn)),
                                            seq_nn,Time)
          colnames(temp_time) <- name_time
          dfTIME              <- rbind(dfTIME,temp_time)
          
          
        }#end for(method in seq_method)
      }#end seq_type_trt
      
      save.image(paste("/Users/santam13/Documents/NYU/ROW/ROW-time-to-event/code/simulations/data_results/simu_n_sample_size.Rdata"))
      
      
      
      ###################################################################################################################################################
      ###################################################################################################################################################
      # Data generation over num covariates
      ###################################################################################################################################################
      ###################################################################################################################################################
      
      beta_ou          <- 0.1 #was  0.5
      beta_tr          <- 0.1 #was 0.5
     
      seq_num_cova    <- c(1,5,10,20,50,100)
      seq_method_bin  <- seq_method_con <- "ROW"
      
      ###################### **************** Over num cova
      rateCens    <- 0.0000000001
      beta_miss   <- 0
      beta_pos    <- 1
      scenario    <- 1 #Continuous treatment with minimal lack of overlap
      for(type_treatment in seq_type_trt){
        for(num_cova in seq_num_cova){
          print(paste(" ###################### **************** Num Covariates ", num_cova, " ****************  ###################### "))
          
          for(i in 1:itera){
            dta <- data_generation_cova(n=n, lambda=lambda, rho=rho, delta = true_haz, 
                                   beta_tr=beta_pos, beta_miss=beta_miss,
                                   beta_ou=beta_ou, rateC=rateCens, 
                                   type_treatment=type_treatment, scenario=scenario,
                                   num_cova=num_cova)
            write.csv(dta, file=paste("~/path",
                                      rateCens,"_type_treatment_",type_treatment,"_lackover_",beta_pos,"_misspe_",beta_miss,"_scenario_",
                                      scenario,"_n_",n,"_","num_cova",num_cova,"_",i,".csv",sep=""), row.names=F)
          }#end for(i in 1:itera)
          
        }#end for(rateCens in seq_rateC)
      }#end seq_type_trt
      
      
      ###################################################################################################################################################
      # Simulations over num covariates
      ###################################################################################################################################################
      ###################################################################################################################################################
      
      seq_method_bin <- seq_method_con <- "ROW"
      
      dfBIAS <- dfMSE <- dfBAL <- dfTIME <- data.frame()
      
      name_bias <- c("Method", "Treatment", "Number Covariates", "Bias")
      name_mse  <- c("Method", "Treatment", "Number Covariates", "MSE")
      name_bal  <- c("Method", "Treatment", "Number Covariates", "Balance")
      name_time <- c("Method", "Treatment", "Number Covariates", "Time")
      
      rateCens    <- 0.0000000001
      beta_miss  <- 0    #No model misspecification
      beta_pos  <- 1 #Moderate lack of overlap for binary treatment
      scenario    <- 1 #Minimal lack of overlap for continuous treatment
      for(type_treatment in seq_type_trt){
        
        if(type_treatment=="binary"){
          seq_method <- seq_method_bin
        }#end if type_treatment
        else{
          seq_method <- seq_method_con
        }#end else
        
        for(method in seq_method){
          print(paste("########################### Method --->", method))
          j <- 0
          for(num_cova in seq_num_cova){
            j <- j+1
            print(paste(" ###################### **************** Num Covariates ", num_cova, " ****************  ###################### "))
            res     <- compute_haz_hat_cova(itera,method,rateCens,type_treatment,beta_pos,beta_miss,n,scenario,delta,num_cova)
            Bias[j] <- res$bias
            MSE[j]  <- res$mse
            Bal[j]  <- mean(res$bal,na.rm=T)
            Time[j] <- res$time
            Pcens[j]<- res$pcens
            
          }#end for(beta_pos in seq_pos){
          
          
          ########################################################################
          #BIAS
          temp_bias           <- data.frame(rep(method,length(seq_num_cova)),
                                            rep(type_treatment,length(seq_num_cova)),
                                            seq_num_cova,Bias)
          colnames(temp_bias) <- name_bias
          dfBIAS              <- rbind(dfBIAS,temp_bias)
          
          ########################################################################
          #MSE
          temp_mse            <- data.frame(rep(method,length(seq_num_cova)),
                                            rep(type_treatment,length(seq_num_cova)),
                                            seq_num_cova,MSE)
          colnames(temp_mse)  <- name_mse
          dfMSE               <- rbind(dfMSE,temp_mse)
          
          ########################################################################
          #BAL
          temp_bal            <- data.frame(rep(method,length(seq_num_cova)),
                                            rep(type_treatment,length(seq_num_cova)),
                                            seq_num_cova,Bal)
          colnames(temp_bal)  <- name_bal
          dfBAL               <- rbind(dfBAL,temp_bal)
          
          ########################################################################
          #TIME
          temp_time           <- data.frame(rep(method,length(seq_num_cova)),
                                            rep(type_treatment,length(seq_num_cova)),
                                            seq_num_cova,Time)
          colnames(temp_time) <- name_time
          dfTIME              <- rbind(dfTIME,temp_time)
          
          
        }#end for(method in seq_method)
      }#end seq_type_trt
      
      save.image(paste("data_results/simu_n",n,"num_cova.Rdata"))
      
      
      
      
      
      
      


      
      
      
