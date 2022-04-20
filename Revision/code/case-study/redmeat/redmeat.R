########################################################################
########################################################################
#***********************************************************************
#
#   In this file, we provide code to evaluate the impact of red meat
#   consumption on time to colon cancer
# 
#***********************************************************************
########################################################################
########################################################################

rm(list=ls())

library(tidyverse)
library(foreign)
library(gurobi)
library(survey)
library(mice)
library(cobalt)
library(SuperLearner)
library(sbw)
library(CBPS)
library(WeightIt)

row <- function(confounders,intervention,delta=0.01){
  
  Xs <- scale(cbind(confounders))
  trs <- as.numeric(scale(intervention))
  n <- dim(confounders)[1]
  
  tol <- 1e-08
  reptimes <- dim(Xs)[2]
  model <- list()
  model$A <- matrix(c(rep(Xs*trs,2),rep(1,n)), nrow = (reptimes*2+1), byrow = T)
  model$rhs <- c(rep(delta,reptimes),rep(-delta,reptimes),1)
  model$modelsense <- "min"
  model$Q <- diag(n)
  model$obj <- rep(1/n,n)
  model$sense <- c(rep("<=",reptimes),rep(">=",reptimes),"=")
  model$lb <- rep(tol, n)
  model$vtypes <- "C"
  params <- list(Presolve = 2, OutputFlag = 1)
  res <- gurobi(model,params)
  return(res)
  
}


data <- read.dta("~/Documents/NYU/ROW/ROW-time-to-event/code/case-study/final_os_12.dta")


data <- as_tibble(data)

data <- data %>% dplyr::select(TIME_Colon, failure_c, totp, f45multi,f45mvmin,ethnic,parity,
                               booph,meno,brca_f2,colon_f2,endo_f2,skin_f2,melan_f2,
                               othca10y,dvt,stroke,mi,diab,hicholrp,osteopor,cvd,cabg,
                               atrialfb,aortican,angina,hip55,smokevr,alcohol,fruits,
                               vegtabls,f60enrgy,syst_bl,dias_bl,bmix_bl,educ,income,
                               time_since_menopause,redmeat)

data <- data %>% mutate(TIME_Colon = TIME_Colon/365)
data <- mutate_all(data, function(x) as.numeric(x))



################################################################################################
#Multiple Imputation

Xm <- data %>% dplyr::select(totp, f45multi,f45mvmin,ethnic,parity,
                             booph,meno,brca_f2,colon_f2,endo_f2,skin_f2,melan_f2,
                             othca10y,dvt,stroke,mi,diab,hicholrp,osteopor,cvd,cabg,
                             atrialfb,aortican,angina,hip55,smokevr,alcohol,fruits,
                             vegtabls,f60enrgy,syst_bl,dias_bl,bmix_bl,educ,income,
                             time_since_menopause,redmeat)

print(sapply(Xm, function(x) sum(is.na(x))/dim(data)[1]))


load("~/Documents/NYU/ROW/ROW-time-to-event/code/case-study/imputed_data_c.Rdata")

data_imputed <- imputed

data_imputed$osteopor <- recode(data_imputed$osteopor, '1'="No", '2'="Yes")
data_imputed$f45multi <- recode(data_imputed$f45multi, '0'="No", '1'="Yes")
data_imputed$f45mvmin <- recode(data_imputed$f45mvmin, '0'="No", '1'="Yes")
data_imputed$brca_f2 <- recode(data_imputed$brca_f2, '0'="No", '1'="Yes")
data_imputed$endo_f2 <- recode(data_imputed$endo_f2, '0'="No", '1'="Yes")
data_imputed$colon_f2 <- recode(data_imputed$colon_f2, '0'="No", '1'="Yes")
data_imputed$skin_f2 <- recode(data_imputed$skin_f2, '0'="No", '1'="Yes")
data_imputed$melan_f2 <- recode(data_imputed$melan_f2, '0'="No", '1'="Yes")
data_imputed$othca10y <- recode(data_imputed$othca10y, '0'="No", '1'="Yes")
data_imputed$dvt <- recode(data_imputed$dvt, '1'="No", '2'="Yes")
data_imputed$stroke <- recode(data_imputed$stroke, '1'="No", '2'="Yes")
data_imputed$mi <- recode(data_imputed$mi, '1'="No", '2'="Yes")
data_imputed$diab <- recode(data_imputed$diab, '0'="No", '1'="Yes")
data_imputed$cvd <- recode(data_imputed$cvd, '1'="No", '2'="Yes")
data_imputed$cabg <- recode(data_imputed$cabg, '1'="No", '2'="Yes")
data_imputed$atrialfb <- recode(data_imputed$atrialfb, '1'="No", '2'="Yes")
data_imputed$aortican <- recode(data_imputed$aortican, '1'="No", '2'="Yes")
data_imputed$angina <- recode(data_imputed$angina, '1'="No", '2'="Yes")
data_imputed$hip55 <- recode(data_imputed$hip55, '1'="No", '2'="Yes")
data_imputed$smokevr <- recode(data_imputed$smokevr, '1'="No", '2'="Yes")
data_imputed$hicholrp <- recode(data_imputed$hicholrp, '1'="No", '2'="Yes")
data_imputed$booph <- recode(data_imputed$booph, '1'="No", '2'="Yes")


data_imputed$TIME_Colon <- data$TIME_Colon
data_imputed$failure_c <- data$failure_c


#Plots

data_imputed_plot <- data_imputed 
confounders_plot <- data_imputed_plot %>% dplyr::select(-c(TIME_Colon, failure_c,redmeat))


#Analysis
data_imputed <- mutate_if(data_imputed,is.character,list(~as.numeric(as.factor(.))))
confounders <- data_imputed %>% dplyr::select(-c(TIME_Colon, failure_c, redmeat))
intervention <- data_imputed$redmeat



## Using alternate variable names



## Using alternate variable names
v1 <- c(totp_Yes = "Hormone Therapy ever",
        f45multi_Yes = "Multivitamine without minerals", 
        f45mvmin_Yes = "Multivitamine with minerals", 
        ethnic = "Ethnicity", 
        parity = "Number of pregnancies", 
        booph_Yes = "Bilateral oophorectomy", 
        meno = "Age at menopause", 
        brca_f2_Yes = "Breast cancer ever", 
        colon_f2_Yes = "Colon cancer ever", 
        endo_f2_Yes = "Endometrial cancer ever", 
        skin_f2_Yes = "Skin cancer ever", 
        melan_f2_Yes = "Melanoma cancer ever",
        othca10y_Yes = "Other cancer past 10 years", 
        dvt_Yes = "Deep vein thrombosis ever", 
        stroke_Yes = "Stroke ever", 
        mi_Yes = "MI ever", 
        diab_Yes = "Diabetes ever", 
        hicholrp_Yes = "High cholesterol pills ever",
        osteopor_Yes = "Osteoporosis ever", 
        cvd_Yes = "Cardiovascular disease ever", 
        cabg_Yes = "Coronary Artery Bypass Graft",
        atrialfb_Yes = "Atrial fibrillation ever", 
        aortican_Yes = "Aortic aneurysm ever", 
        angina_Yes ="Angina", 
        hip55_Yes ="Hip fracture age 55 or older",
        smokevr_Yes ="Smoked >100 cigarettes ever", 
        alcohol ="Alcohol intake", 
        fruits = "Fruits, med serv/day",
        vegtabls = "Vegetables, med serv/day", 
        f60enrgy = "Dietary Energy (kcal)", 
        syst_bl = "Systolic blood pressure",
        dias_bl ="Diastolic blood pressure", 
        bmix_bl = "Body Mass Index", 
        educ = "Education", 
        income = "Income",
        time_since_menopause = "Time since menopause")



####################################################################################################
####################################################################################################
####################################################################################################
#ROW

timeROW    <- tryCatch(  system.time(row_w <- row(confounders,intervention, delta = 0.001))
                         , error=function(e) NULL)#endtryCatch

n <- dim(confounders)[1]
summary(row_w$x*n)

lprow <- love.plot(bal.tab(intervention~confounders_plot,
                           weights=row_w$x,
                           method="weighting",
                           s.d.denom="all"),
                   colors = c("grey", "black"),
                   var.order = "unadjusted",
                   shapes=c("circle","square"),
                   stars = "raw", var.names = v1, legend="none", position=1, abs=TRUE)

lprow <- lprow + theme(axis.text.x = element_text(size = 14),
                         axis.text.y = element_text(size = 12, angle = 0),
                         text = element_text(size = 14))

am_row <- bal.tab(intervention~confounders,
                  weights=row_w$x,
                  method="weighting")$Balance$Corr.Adj

surveyDesign  <- svydesign(ids=1:n,weights=~row_w$x,data=data_imputed)
modelw   <- svycoxph(Surv(TIME_Colon, failure_c) ~ redmeat,design=surveyDesign) 
summary(modelw)



############
#Test for Reviewer, let's remove the most imbalanced covariates such as energy and bmi from the analysis


confounders_small <- confounders %>% dplyr::select(-c(f60enrgy,bmix_bl,educ,fruits))

timeROW    <- tryCatch(  system.time(row_w <- row(confounders_small,intervention, delta = 0.001))
                         , error=function(e) NULL)#endtryCatch

n <- dim(confounders)[1]
summary(row_w$x*n)


surveyDesign  <- svydesign(ids=1:n,weights=~row_w$x,data=data_imputed)
modelw_small   <- svycoxph(Surv(TIME_Colon, failure_c) ~ redmeat,design=surveyDesign) 
summary(modelw_small)


####################################################################################################
####################################################################################################
####################################################################################################
#CBPS

timeCBPS    <- tryCatch(  system.time(cbps_w  <- CBPS(scale(intervention) ~ scale(totp)+scale(f45multi)+scale(f45mvmin)+scale(ethnic)+scale(parity)+scale(booph)+scale(meno)+scale(brca_f2)+scale(colon_f2)+scale(endo_f2)+scale(skin_f2)+scale(melan_f2)+scale(
  othca10y)+scale(dvt)+scale(stroke)+scale(mi)+scale(diab)+scale(hicholrp)+scale(osteopor)+scale(cvd)+
    scale(cabg)+scale(atrialfb)+scale(aortican)+scale(angina)+scale(hip55)+scale(smokevr)+scale(alcohol)+
    scale(fruits)+scale(vegtabls)+scale(f60enrgy)+scale(syst_bl)+scale(dias_bl)+scale(bmix_bl)+scale(educ)+
    scale(income)+scale(time_since_menopause),method = "exact", data=data_imputed, ATT=0)$weights)
  , error=function(e) NULL)#endtryCatch

lpcbps <- love.plot(bal.tab(intervention~confounders_plot,
                            weights=cbps_w,
                            method="weighting",
                            s.d.denom="all"),
                    colors = c("grey", "black"),
                    var.order = "unadjusted",
                    shapes=c("circle","square"),
                    stars = "raw", var.names = v, legend="none", position=1, abs=TRUE)


am_cbps <- bal.tab(intervention~confounders,
                   weights=cbps_w,
                   method="weighting")$Balance$Corr.Adj

surveyDesign  <- svydesign(ids=1:n,weights=~cbps_w,data=data_imputed)
modelcbps   <- svycoxph(Surv(TIME_Colon, failure_c) ~ redmeat,design=surveyDesign) 
summary(modelcbps)


####################################################################################################
####################################################################################################
####################################################################################################
#IPW


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
  
  
  m1 <- median(X$meno)
  m2 <- median(X$fruits)
  m3 <- median(X$vegtabls)
  m4 <- median(X$f60enrgy)
  m5 <- median(X$syst_bl)
  m6 <- median(X$dias_bl)
  m7 <- median(X$bmix_bl)
  m8 <- median(X$time_since_menopause)
  
  
  fit.glm <- glm(Y ~ bs(meno,knots=m1) + 
                   bs(fruits,knots=m2) + 
                   bs(vegtabls,knots=m3) +
                   bs(f60enrgy,knots=m4) + 
                   bs(syst_bl,knots=m5) + 
                   bs(dias_bl,knots=m6) +
                   bs(bmix_bl,knots=m7) +
                   bs(time_since_menopause,knots=m8) 
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
  
  
  
  fit.glm <- glm(Y ~ fp(meno,df=2) + 
                   fp(fruits,df=2) + 
                   fp(vegtabls,df=2) +
                   fp(f60enrgy,df=2) + 
                   fp(syst_bl,df=2) + 
                   fp(dias_bl,df=2) +
                   fp(bmix_bl,df=2) +
                   fp(time_since_menopause,df=2) 
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

sl.lib_gps <- c('SL.glm',"SL.ranger",'SL.glm.FP','SL.glm.BS')

timeipw    <- tryCatch(  system.time(gps <- weightit(intervention ~ confounders , data = data_imputed,
                                                     method = "super", density = "dt_2",stabilize = FALSE,
                                                     SL.library = sl.lib_gps, family="gaussian"))
                         , error=function(e) NULL)#endtryCatch

ipw <- gps$weights



lpipw1 <- love.plot(bal.tab(intervention~confounders_plot,
                            weights=ipw,
                            method="weighting",
                            s.d.denom="all"),
                    colors = c("grey", "black"),
                    var.order = "unadjusted",
                    shapes=c("circle","square"),
                    stars = "raw", var.names = v1, legend="none", position=1, abs=TRUE)

am_ipw <- bal.tab(intervention~confounders,
                  weights=ipw,
                  method="weighting")$Balance$Corr.Adj

surveyDesign  <- svydesign(ids=1:n,weights=~ipw,data=data_imputed)
modelipw   <- svycoxph(Surv(TIME_Colon, failure_c) ~ redmeat,design=surveyDesign) 
summary(modelipw)

summary(coxph(Surv(TIME_Colon, failure_c) ~ redmeat,weights=ipw,data=data_imputed))



####################################################################################################
####################################################################################################
####################################################################################################
#IPW-trunc

ipw_T <- ipw

QQ <- 0.99

ipw_T[which(ipw_T>=quantile(ipw_T,QQ))] <- quantile(gps$weights,c(QQ))
ipw_T[which(ipw_T<=quantile(ipw_T,1-QQ))] <- quantile(gps$weights,c(1-QQ))


lpipw1_T <- love.plot(bal.tab(intervention~confounders_plot,
                            weights=ipw_T,
                            method="weighting",
                            s.d.denom="all"),
                    colors = c("grey", "black"),
                    var.order = "unadjusted",
                    shapes=c("circle","square"),
                    stars = "raw", var.names = v, legend="none", position=1, abs=TRUE)

am_ipw_T <- bal.tab(intervention~confounders,
                  weights=ipw_T,
                  method="weighting")$Balance$Corr.Adj

surveyDesign  <- svydesign(ids=1:n,weights=~ipw_T,data=data_imputed)
modelipw_T   <- svycoxph(Surv(TIME_Colon, failure_c) ~ redmeat,design=surveyDesign) 
summary(modelipw_T)

summary(coxph(Surv(TIME_Colon, failure_c) ~ redmeat,weights=ipw_T,data=data_imputed))


####################################################################################################
####################################################################################################
####################################################################################################
#IPW-stab

stabilize <- TRUE
truncate <- FALSE
qt <- 0.1


sl.lib_gps <- c('SL.glm',"SL.ranger",'SL.glm.FP','SL.glm.BS')

timeipw_S    <- tryCatch(  system.time(gps_S <- weightit(intervention ~ confounders , data = data_imputed,
                                                     method = "super", density = "dt_2", stabilize = stabilize,
                                                     SL.library = sl.lib_gps, family="gaussian"))
                         , error=function(e) NULL)#endtryCatch

ipw_S <- gps_S$weights


lpipw1 <- love.plot(bal.tab(intervention~confounders_plot,
                            weights=ipw_S,
                            method="weighting",
                            s.d.denom="all"),
                    colors = c("grey", "black"),
                    var.order = "unadjusted",
                    shapes=c("circle","square"),
                    stars = "raw", var.names = v, legend="none", position=1, abs=TRUE)

am_ipw_S <- bal.tab(intervention~confounders,
                  weights=ipw_S,
                  method="weighting")$Balance$Corr.Adj

surveyDesign  <- svydesign(ids=1:n,weights=~ipw_S,data=data_imputed)
modelipw_S   <- svycoxph(Surv(TIME_Colon, failure_c) ~ redmeat,design=surveyDesign) 
summary(modelipw_S)


summary(coxph(Surv(TIME_Colon, failure_c) ~ redmeat,weights=ipw_S,data=data_imputed))


####################################################################################################
####################################################################################################
####################################################################################################
#BalSL

sl.lib_gps <- c('SL.glm',"SL.ranger")

timeBalSL    <- tryCatch(  system.time(ps_BalSL <- weightit(intervention ~ confounders , data = data_imputed,
                                                            method = "super", density = "dt_2",
                                                            SL.library = sl.lib_gps, family="gaussian",
                                                            SL.method = "method.balance",
                                                            stop.method = "p.mean"))
                           , error=function(e) NULL)#endtryCatch


lpBalSL <- love.plot(bal.tab(intervention~confounders_plot,
                             weights=ps_BalSL$weights,
                             method="weighting",
                             s.d.denom="all"),
                     colors = c("grey", "black"),
                     var.order = "unadjusted",
                     shapes=c("circle","square"),
                     stars = "raw", var.names = v, legend="none", position=1, abs=TRUE)

am_bsl <- bal.tab(intervention~confounders,
                  weights=ps_BalSL$weights,
                  method="weighting")$Balance$Corr.Adj

surveyDesign  <- svydesign(ids=1:n,weights=~ps_BalSL$weights,data=data_imputed)
modelbalSL   <- svycoxph(Surv(TIME_Colon, failure_c) ~ redmeat,design=surveyDesign) 
summary(modelbalSL)




####################################################################################################
####################################################################################################
####################################################################################################
#GBM

timegbm    <- tryCatch(  system.time(ps_gbm <- weightit(intervention ~ confounders , data = data_imputed,
                                                        method = "gbm", density = "dnorm",
                                                        family="gaussian",
                                                        stop.method = "p.mean"))
                         , error=function(e) NULL)#endtryCatch


lpgbm <- love.plot(bal.tab(intervention~confounders_plot,
                           weights=ps_gbm$weights,
                           method="weighting",
                           s.d.denom="all"),
                   colors = c("grey", "black"),
                   var.order = "unadjusted",
                   shapes=c("circle","square"),
                   stars = "raw", var.names = v, legend="none", position=1, abs=TRUE)

am_gbm <- bal.tab(intervention~confounders,
                  weights=ps_gbm$weights,
                  method="weighting")$Balance$Corr.Adj

surveyDesign  <- svydesign(ids=1:n,weights=~ps_gbm$weights,data=data_imputed)
modelgbm   <- svycoxph(Surv(TIME_Colon, failure_c) ~ redmeat,design=surveyDesign) 
summary(modelgbm)

####################################################################################################
####################################################################################################
####################################################################################################
#npCBPS

timenpcbps    <- tryCatch(  system.time(npcbps_w  <- npCBPS(scale(intervention) ~ scale(totp)+scale(f45multi)+scale(f45mvmin)+scale(ethnic)+scale(parity)+scale(booph)+scale(meno)+scale(brca_f2)+scale(colon_f2)+scale(endo_f2)+scale(skin_f2)+scale(melan_f2)+scale(
  othca10y)+scale(dvt)+scale(stroke)+scale(mi)+scale(diab)+scale(hicholrp)+scale(osteopor)+scale(cvd)+
    scale(cabg)+scale(atrialfb)+scale(aortican)+scale(angina)+scale(hip55)+scale(smokevr)+scale(alcohol)+
    scale(fruits)+scale(vegtabls)+scale(f60enrgy)+scale(syst_bl)+scale(dias_bl)+scale(bmix_bl)+scale(educ)+
    scale(income)+scale(time_since_menopause),method = "exact", data=data_imputed, ATT=0)$weights )
  , error=function(e) NULL)#endtryCatch

lpnpcbps <- love.plot(bal.tab(intervention~confounders_plot,
                              weights=npcbps_w,
                              method="weighting",
                              s.d.denom="all"),
                      colors = c("grey", "black"),
                      var.order = "unadjusted",
                      shapes=c("circle","square"),
                      stars = "raw", var.names = v)

am_npcbps <- bal.tab(intervention~confounders,
                     weights=npcbps_w,
                     method="weighting")$Balance$Corr.Adj

surveyDesign  <- svydesign(ids=1:n,weights=~npcbps_w,data=data_imputed)
modelnpcbps   <- svycoxph(Surv(TIME_Colon, failure_c) ~ redmeat,design=surveyDesign) 
summary(modelnpcbps)


####################################################################################################
####################################################################################################
####################################################################################################
#Outcome modelling
timeOM    <- tryCatch(  system.time(model_om <- coxph(Surv(TIME_Colon, failure_c)~intervention + totp+f45multi+f45mvmin+ethnic+parity+booph+meno+
                                                        brca_f2+colon_f2+endo_f2+skin_f2+melan_f2+othca10y+dvt+stroke+mi+diab+hicholrp+osteopor+cvd+
                                                        cabg+atrialfb+aortican+angina+hip55+smokevr+alcohol+
                                                        fruits+vegtabls+f60enrgy+syst_bl+dias_bl+bmix_bl+educ+
                                                        income+time_since_menopause,data=data_imputed))
                        , error=function(e) NULL)#endtryCatch
exp(summary(model_om)$coeff[1])


####################################################################################################
####################################################################################################
####################################################################################################
#Outcome modelling
timeNA    <- tryCatch(  system.time(model_na <- coxph(Surv(TIME_Colon, failure_c)~intervention,data=data_imputed))
                        , error=function(e) NULL)#endtryCatch
exp(summary(model_na)$coeff[1])


am_una <- bal.tab(intervention~confounders)$Balance$Corr.Un

save.image("case_study_c2.Rdata") 



load("case_study_c2.Rdata")


table_1 <- rbind(c(round(exp(modelw$coefficients[1]),2),round(exp(confint(modelw)),2), round(timeROW[1]+timeROW[2],1), 
                   round(median(abs(am_row)),4), round(min(abs(am_row)),4), round(max(abs(am_row)),4) ),
                 c(round(exp(model_na$coefficients[1]),2),round(exp(confint(model_na)),2), round(timeNA[1]+timeNA[2],1), 
                   c(0,0,0)),
                 c(round(exp(modelbalSL$coefficients[1]),2),round(exp(confint(modelbalSL)),2), round(timeBalSL[1]+timeBalSL[2],1) , 
                   round(median(abs(am_bsl)),4), round(min(abs(am_bsl)),4), round(max(abs(am_bsl)),4) ),
                 c(round(exp(modelgbm$coefficients[1]),2),round(exp(confint(modelgbm)),2),round(timegbm[1]+timegbm[2],1), 
                   round(median(abs(am_gbm)),4), round(min(abs(am_gbm)),4), round(max(abs(am_gbm)),4) ),
                 c(round(exp(model_om$coefficients[1]),2),round(exp(confint(model_om))[1,],2),round(timeOM[1]+timeOM[2],1), 
                   c(0,0,0)),
                 c(round(exp(modelipw$coefficients[1]),2),round(exp(confint(modelipw)),2),round(timeipw[1]+timeipw[2],1), 
                   round(median(abs(am_ipw)),4), round(min(abs(am_ipw)),4), round(max(abs(am_ipw)),4) ),
                 c(round(exp(modelcbps$coefficients[1]),2),round(exp(confint(modelcbps)),2),round(timeCBPS[1]+timeCBPS[2],1), 
                   round(median(abs(am_cbps)),4), round(min(abs(am_cbps)),4), round(max(abs(am_cbps)),4) ),
                 c(round(exp(modelnpcbps$coefficients[1]),2),round(exp(confint(modelnpcbps)),2),round(timenpcbps[1]+timenpcbps[2],1), 
                   round(median(abs(am_npcbps)),4), round(min(abs(am_npcbps)),4), round(max(abs(am_npcbps)),4) )
                   )


rownames(table_1) <- c("ROW","Naive","BalSL","GBM","OM","IPW","CBPS","npCBPS")
colnames(table_1) <- c("HR","2.5%","97.5%","Time in seconds", "Median absolute balance", "Min absolute balance", "Max absolute balance")


table_1

write.csv(table_1,'table1c.csv')


#########################################################################################################
# Extract the legend. Returns a gtable
grid.newpage()
pushViewport(viewport(layout = grid.layout(1, 1, heights = unit(c(0.5, 5), "null"))))  


print(lprow, vp = viewport(layout.pos.row = 1, layout.pos.col = 1))







