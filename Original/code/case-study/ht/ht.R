########################################################################
########################################################################
#***********************************************************************
#
#   In this file, we provide code to evaluate the impact of hormone therapy
#   estrogin plus progestin on time to CHD
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
library(Matching)
library(ggpubr)
library(grid)
library(survminer)

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


data <- read.dta("~/path/final_os_12.dta")

data <- as_tibble(data)

data <- data %>% dplyr::select(TIME_CHD, failure, totp, f45multi,f45mvmin,ethnic,parity,
                               booph,meno,brca_f2,colon_f2,endo_f2,skin_f2,melan_f2,
                               othca10y,dvt,stroke,mi,diab,hicholrp,osteopor,cvd,cabg,
                               atrialfb,aortican,angina,hip55,smokevr,alcohol,fruits,
                               vegtabls,f60enrgy,syst_bl,dias_bl,bmix_bl,educ,income,
                               time_since_menopause)

data <- data %>% mutate(TIME_CHD = TIME_CHD/365)
data <- mutate_all(data, function(x) as.numeric(x))


################################################################################################
#Multiple Imputation

Xm <- data %>% dplyr::select(totp, f45multi,f45mvmin,ethnic,parity,
                             booph,meno,brca_f2,colon_f2,endo_f2,skin_f2,melan_f2,
                             othca10y,dvt,stroke,mi,diab,hicholrp,osteopor,cvd,cabg,
                             atrialfb,aortican,angina,hip55,smokevr,alcohol,fruits,
                             vegtabls,f60enrgy,syst_bl,dias_bl,bmix_bl,educ,income,
                             time_since_menopause)

print(sapply(Xm, function(x) sum(is.na(x))/dim(data)[1]))

init  <-  mice(Xm, maxit=0)
meth  <-  init$method
predM <-  init$predictorMatrix

imputed <- mice(Xm, method="pmm", predictorMatrix=predM, m=5)
imputed <- mice::complete(imputed)

save.image("imputed_data.Rdata")
load("imputed_data.Rdata")


data_imputed <- imputed

data_imputed$osteopor <- recode(data_imputed$osteopor, '1'="No", '2'="Yes")
data_imputed$f45multi <- recode(data_imputed$f45multi, '1'="No", '2'="Yes")
data_imputed$f45mvmin <- recode(data_imputed$f45mvmin, '1'="No", '2'="Yes")
data_imputed$brca_f2 <- recode(data_imputed$brca_f2, '1'="No", '2'="Yes")
data_imputed$endo_f2 <- recode(data_imputed$endo_f2, '1'="No", '2'="Yes")
data_imputed$colon_f2 <- recode(data_imputed$colon_f2, '1'="No", '2'="Yes")
data_imputed$skin_f2 <- recode(data_imputed$skin_f2, '1'="No", '2'="Yes")
data_imputed$melan_f2 <- recode(data_imputed$melan_f2, '1'="No", '2'="Yes")
data_imputed$othca10y <- recode(data_imputed$othca10y, '1'="No", '2'="Yes")
data_imputed$dvt <- recode(data_imputed$dvt, '1'="No", '2'="Yes")
data_imputed$stroke <- recode(data_imputed$stroke, '1'="No", '2'="Yes")
data_imputed$mi <- recode(data_imputed$mi, '1'="No", '2'="Yes")
data_imputed$diab <- recode(data_imputed$diab, '1'="No", '2'="Yes")
data_imputed$cvd <- recode(data_imputed$cvd, '1'="No", '2'="Yes")
data_imputed$cabg <- recode(data_imputed$cabg, '1'="No", '2'="Yes")
data_imputed$atrialfb <- recode(data_imputed$atrialfb, '1'="No", '2'="Yes")
data_imputed$aortican <- recode(data_imputed$aortican, '1'="No", '2'="Yes")
data_imputed$angina <- recode(data_imputed$angina, '1'="No", '2'="Yes")
data_imputed$hip55 <- recode(data_imputed$hip55, '1'="No", '2'="Yes")
data_imputed$smokevr <- recode(data_imputed$smokevr, '1'="No", '2'="Yes")
data_imputed$hicholrp <- recode(data_imputed$hicholrp, '1'="No", '2'="Yes")
data_imputed$booph <- recode(data_imputed$booph, '1'="No", '2'="Yes")


data_imputed$TIME_CHD <- data$TIME_CHD
data_imputed$failure <- data$failure
data_imputed$time_since_menopause <- data$time_since_menopause


#Plots

data_imputed_1_plot <- data_imputed %>% filter(time_since_menopause>0 & time_since_menopause<=10)
confounders_1_plot <- data_imputed_1_plot %>% dplyr::select(-c(TIME_CHD, failure,totp,time_since_menopause))


data_imputed_2_plot <- data_imputed %>% filter(time_since_menopause>10 & time_since_menopause<=20)
confounders_2_plot <- data_imputed_2_plot %>% dplyr::select(-c(TIME_CHD, failure,totp, time_since_menopause))


data_imputed_3_plot <- data_imputed %>% filter(time_since_menopause>20)
confounders_3_plot <- data_imputed_3_plot %>% dplyr::select(-c(TIME_CHD, failure,totp,time_since_menopause))



#Analysis
data_imputed <- mutate_if(data_imputed,is.character,list(~as.numeric(as.factor(.))))

data_imputed_1 <- data_imputed %>% filter(time_since_menopause>0 & time_since_menopause<=10)
confounders_1 <- data_imputed_1 %>% dplyr::select(-c(TIME_CHD, failure,totp, time_since_menopause))
intervention_1 <- data_imputed_1$totp

data_imputed_2 <- data_imputed %>% filter(time_since_menopause>10 & time_since_menopause<=20)
confounders_2 <- data_imputed_2 %>% dplyr::select(-c(TIME_CHD, failure,totp, time_since_menopause))
intervention_2 <- data_imputed_2$totp

data_imputed_3 <- data_imputed %>% filter(time_since_menopause>20)
confounders_3 <- data_imputed_3 %>% dplyr::select(-c(TIME_CHD, failure,totp,time_since_menopause))
intervention_3 <- data_imputed_3$totp



## Using alternate variable names
v1 <- c(f45multi_Yes = "Multivitamine without minerals", 
        f45mvmin_Yes = "Multivitamine with minerals", 
        ethnic = "Ethnicity", 
        parity = "Nuber of pregnancies", 
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
        income = "Income")



####################################################################################################
####################################################################################################
####################################################################################################
#ROW

###########################################################################
#time_since_menopause 0-10
timeROW1    <- tryCatch(  system.time(row_w <- row(confounders_1,intervention_1, delta = 0.0001))
                          , error=function(e) NULL)#endtryCatch

n <- dim(confounders_1)[1]
summary(row_w$x*n)

am_row_1 <- bal.tab(intervention_1~confounders_1,
            weights=row_w$x,
            method="weighting",
            s.d.denom="pooled")$Balance$Diff.Adj


lprow1 <- love.plot(bal.tab(intervention_1~confounders_1_plot,
                            weights=row_w$x,
                            method="weighting",
                            s.d.denom="pooled"),
                    colors = c("grey", "black"),
                    var.order = "unadjusted",
                    shapes=c("circle","square"),
                    stars = "none", var.names = v1, 
                    legend="none", position=1, abs=TRUE,
                    sample.names = c("Unweighted", "Weighted")
                    )

lprow1 <- lprow1 + theme(axis.text.x = element_text(size = 14),
                         axis.text.y = element_text(size = 12, angle = 0),
               text = element_text(size = 14))

surveyDesign_1  <- svydesign(ids=1:n,weights=~row_w$x,data=data_imputed_1)
modelw_1   <- svycoxph(Surv(TIME_CHD, failure) ~ totp,design=surveyDesign_1) 
summary(modelw_1)


cox.zph(modelw_1)

km1 <- ggsurvplot(survfit(Surv(TIME_CHD/365, failure) ~ totp,data=data_imputed_1,weights=(row_w$x*n)),
                  risk.table=TRUE, pval = FALSE, pval.method=FALSE, pval.coord= c(100,0.75),
                  pval.method.coord= c(100,0.80), legend.title = "HT",legend.labs = c("No", "Yes"),
                  cumevents=FALSE,cumcensor=FALSE, ylim= c(0.5,1), xlab="Time in years since baseline",              palette = c("black", "grey"),
                  ggtheme = theme_bw(), legend="right")

km1$plot  <- km1$plot + theme(axis.text.x = element_text(size = 14),
                              axis.text.y = element_text(size = 12, angle = 0),
                              text = element_text(size = 14))


###########################################################################
#time_since_menopause 11-20
timeROW2    <- tryCatch(  system.time(row_w <- row(confounders_2,intervention_2, delta = 0.0001))
                          , error=function(e) NULL)#endtryCatch

n <- dim(confounders_2)[1]
summary(row_w$x*n)

am_row_2 <- bal.tab(intervention_2~confounders_2_plot,
            weights=row_w$x,
            method="weighting",
            s.d.denom="pooled")$Balance$Diff.Adj

lprow2 <- love.plot(bal.tab(intervention_2~confounders_2_plot,
                            weights=row_w$x,
                            method="weighting",
                            s.d.denom="pooled"),
                    colors = c("grey", "black"),
                    var.order = "unadjusted",
                    shapes=c("circle","square"),
                    stars = "none", var.names = v1, legend="none", position=1, abs=TRUE)

lprow2 <- lprow2 + theme(axis.text = element_text(size = 12),
                         text = element_text(size = 14))

surveyDesign_2  <- svydesign(ids=1:n,weights=~row_w$x,data=data_imputed_2)
modelw_2   <- svycoxph(Surv(TIME_CHD, failure) ~ totp,design=surveyDesign_2) 
summary(modelw_2)

cox.zph(modelw_2)

km2 <- ggsurvplot(survfit(Surv(TIME_CHD/365, failure) ~ totp,data=data_imputed_2,weights=(row_w$x*n)),
                  risk.table=TRUE, pval = FALSE, pval.method=FALSE, pval.coord= c(100,0.75),
                  pval.method.coord= c(100,0.80), legend.title = "HT",legend.labs = c("No", "Yes"),
                  cumevents=FALSE,cumcensor=FALSE, ylim= c(0.5,1), xlab="Time in years since baseline",              palette = c("black", "grey"),
                  ggtheme = theme_bw(), legend="right")

km2$plot <- km2$plot + theme(axis.text.x = element_text(size = 14),
                             axis.text.y = element_text(size = 12, angle = 0),
                             text = element_text(size = 14))

###########################################################################
#time_since_menopause 20+
timeROW3    <- tryCatch(  system.time(row_w <- row(confounders_3,intervention_3, delta = 0.0001))
                          , error=function(e) NULL)#endtryCatch

n <- dim(confounders_3)[1]
summary(row_w$x*n)

am_row_3 <- bal.tab(intervention_3~confounders_3_plot,
            weights=row_w$x,
            method="weighting",
            s.d.denom="pooled")$Balance$Diff.Adj

lprow3 <- love.plot(bal.tab(intervention_3~confounders_3_plot,
                            weights=row_w$x,
                            method="weighting",
                            s.d.denom="pooled"),
                    colors = c("grey", "black"),
                    var.order = "unadjusted",
                    shapes=c("circle","square"),
                    stars = "none", var.names = v1, legend="none", position=1, abs=TRUE)

lprow3 <- lprow3 + theme(axis.text = element_text(size = 12),
                         text = element_text(size = 14))

surveyDesign_3  <- svydesign(ids=1:n,weights=~row_w$x,data=data_imputed_3)
modelw_3   <- svycoxph(Surv(TIME_CHD, failure) ~ totp,design=surveyDesign_3) 
summary(modelw_3)

cox.zph(modelw_3)

km3 <- ggsurvplot(survfit(Surv(TIME_CHD/365, failure) ~ totp,data=data_imputed_3,weights=(row_w$x*n)),
                  risk.table=TRUE, pval = FALSE, pval.method=FALSE, pval.coord= c(100,0.75),
                  pval.method.coord= c(100,0.80), legend.title = "HT",legend.labs = c("No", "Yes"),
                  cumevents=FALSE,cumcensor=FALSE, ylim= c(0.5,1), xlab="Time in years since baseline",
                  palette = c("black", "grey"),
                  ggtheme = theme_bw(), legend="right")

km3$plot <- km3$plot + theme(axis.text.x = element_text(size = 14),
                             axis.text.y = element_text(size = 12, angle = 0),
                             text = element_text(size = 14))

####################################################################################################
####################################################################################################
####################################################################################################
#Naive
timeNa1    <- tryCatch(  system.time(modeln_1 <- coxph(Surv(TIME_CHD,failure)~totp,data=data_imputed_1))
                         , error=function(e) NULL)#endtryCatch
summary(modeln_1)
timeNa2    <- tryCatch(  system.time(modeln_2 <- coxph(Surv(TIME_CHD,failure)~totp,data=data_imputed_2))
                         , error=function(e) NULL)#endtryCatch
summary(modeln_2)
timeNa3    <- tryCatch(  system.time(modeln_3 <- coxph(Surv(TIME_CHD,failure)~totp,data=data_imputed_3))
                         , error=function(e) NULL)#endtryCatch
summary(modeln_3)


####################################################################################################
####################################################################################################
####################################################################################################
#Outcome modelling
timeOM1    <- tryCatch(  system.time(modelm_1 <- coxph(Surv(TIME_CHD,failure)~totp+ f45multi+f45mvmin+ethnic+parity+
                                                         booph+meno+brca_f2+colon_f2+endo_f2+skin_f2+melan_f2+
                                                         othca10y+dvt+stroke+mi+diab+hicholrp+osteopor+cvd+cabg+
                                                         atrialfb+aortican+angina+hip55+smokevr+alcohol+fruits+
                                                         vegtabls+f60enrgy+syst_bl+dias_bl+bmix_bl+educ+income,data=data_imputed_1))
                         , error=function(e) NULL)#endtryCatch
exp(summary(modelm_1)$coeff[1])

timeOM2    <- tryCatch(  system.time(modelm_2 <- coxph(Surv(TIME_CHD,failure)~totp+ f45multi+f45mvmin+ethnic+parity+
                                                         booph+meno+brca_f2+colon_f2+endo_f2+skin_f2+melan_f2+
                                                         othca10y+dvt+stroke+mi+diab+hicholrp+osteopor+cvd+cabg+
                                                         atrialfb+aortican+angina+hip55+smokevr+alcohol+fruits+
                                                         vegtabls+f60enrgy+syst_bl+dias_bl+bmix_bl+educ+income,data=data_imputed_2))
                         , error=function(e) NULL)#endtryCatch
exp(summary(modelm_2)$coeff[1])


timeOM3    <- tryCatch(  system.time(modelm_3 <- coxph(Surv(TIME_CHD,failure)~totp+f45multi+f45mvmin+ethnic+parity+
                                                         booph+meno+brca_f2+colon_f2+endo_f2+skin_f2+melan_f2+
                                                         othca10y+dvt+stroke+mi+diab+hicholrp+osteopor+cvd+cabg+
                                                         atrialfb+aortican+angina+hip55+smokevr+alcohol+fruits+
                                                         vegtabls+f60enrgy+syst_bl+dias_bl+bmix_bl+educ+income,data=data_imputed_3))
                         , error=function(e) NULL)#endtryCatch
exp(summary(modelm_3)$coeff[1])


####################################################################################################
####################################################################################################
####################################################################################################
#IPW

sl.lib_ps  <- c('SL.glm','SL.glmnet',"SL.randomForest")

#Balancing covariates between treatment groups (binary)
Tr    <- abs(intervention_1-1)
timeIPW1    <- tryCatch(  system.time( SL    <- SuperLearner(Tr, confounders_1 ,
                                                             SL.library = sl.lib_ps, 
                                                             family = "binomial"))
                          , error=function(e) NULL)#endtryCatch
ps    <- SL$SL.predict
ipw   <- ( Tr/ps + (1-Tr)/(1-ps)  )
summary(ipw)

n <- dim(confounders_1)[1]
am_ipw_1 <- bal.tab(intervention_1~confounders_1,
            weights=ipw,
            method="weighting",
            s.d.denom="pooled")$Balance$Diff.Adj

lpipw1 <- love.plot(bal.tab(intervention_1~confounders_1_plot,
                            weights=ipw,
                            method="weighting",
                            s.d.denom="pooled"),
                    colors = c("grey", "black"),
                    var.order = "unadjusted",
                    shapes=c("circle","square"),
                    stars = "raw", var.names = v, legend="none", position=1, abs=TRUE)

surveyDesignipw_1  <- svydesign(ids=1:n,weights=~ipw,data=data_imputed_1)
modelipw_1   <- svycoxph(Surv(TIME_CHD, failure) ~ totp,design=surveyDesignipw_1) 
summary(modelipw_1)



#Balancing covariates between treatment groups (binary)
Tr    <- abs(intervention_2-1)
timeIPW2    <- tryCatch(  system.time( SL    <- SuperLearner(Tr, confounders_2 ,
                                                             SL.library = sl.lib_ps, 
                                                             family = "binomial"))
                          , error=function(e) NULL)#endtryCatch
ps    <- SL$SL.predict
ipw   <- ( Tr/ps + (1-Tr)/(1-ps)  )
summary(ipw)

n <- dim(confounders_2)[1]
am_ipw_2 <- bal.tab(intervention_2~confounders_2,
            weights=ipw,
            method="weighting",
            s.d.denom="pooled")$Balance$Diff.Adj

lpipw2 <- love.plot(bal.tab(intervention_2~confounders_2_plot,
                            weights=ipw,
                            method="weighting",
                            s.d.denom="pooled"),
                    colors = c("grey", "black"),
                    var.order = "unadjusted",
                    shapes=c("circle","square"),
                    stars = "raw", var.names = v, legend="none", position=1, abs=TRUE)

surveyDesignipw_2  <- svydesign(ids=1:n,weights=~ipw,data=data_imputed_2)
modelipw_2   <- svycoxph(Surv(TIME_CHD, failure) ~ totp,design=surveyDesignipw_2) 
summary(modelipw_2)



#Balancing covariates between treatment groups (binary)
Tr    <- abs(intervention_3-1)
timeIPW3    <- tryCatch(  system.time( SL    <- SuperLearner(Tr, confounders_3 ,
                                                             SL.library = sl.lib_ps, 
                                                             family = "binomial"))
                          , error=function(e) NULL)#endtryCatch
ps    <- SL$SL.predict
ipw   <- ( Tr/ps + (1-Tr)/(1-ps)  )
summary(ipw)

n <- dim(confounders_3)[1]
am_ipw_3 <- bal.tab(intervention_3~confounders_3_plot,
            weights=ipw,
            method="weighting",
            s.d.denom="pooled")$Balance$Diff.Adj

lpipw3 <- love.plot(bal.tab(intervention_3~confounders_3,
                            weights=ipw,
                            method="weighting",
                            s.d.denom="pooled"),
                    colors = c("grey", "black"),
                    var.order = "unadjusted",
                    shapes=c("circle","square"),
                    stars = "raw", var.names = v, legend="none", position=1, abs=TRUE)

surveyDesignipw_3  <- svydesign(ids=1:n,weights=~ipw,data=data_imputed_3)
modelipw_3   <- svycoxph(Surv(TIME_CHD, failure) ~ totp,design=surveyDesignipw_3) 
summary(modelipw_3)



####################################################################################################
####################################################################################################
####################################################################################################
#CBPS
timeCBPS1    <- tryCatch(  system.time( cbps_w  <- CBPS(as.factor(intervention_1) ~ scale(f45multi)+scale(f45mvmin)+scale(ethnic)+scale(parity)+scale(booph)+scale(meno)+scale(brca_f2)+scale(colon_f2)+scale(endo_f2)+scale(skin_f2)+scale(melan_f2)+scale(
  othca10y)+scale(dvt)+scale(stroke)+scale(mi)+scale(diab)+scale(hicholrp)+scale(osteopor)+scale(cvd)+
    scale(cabg)+scale(atrialfb)+scale(aortican)+scale(angina)+scale(hip55)+scale(smokevr)+scale(alcohol)+
    scale(fruits)+scale(vegtabls)+scale(f60enrgy)+scale(syst_bl)+scale(dias_bl)+scale(bmix_bl)+scale(educ)+
    scale(income),method = "exact", data=data_imputed_1, ATT=0)$weights )
  , error=function(e) NULL)#endtryCatch


n <- dim(confounders_1)[1]
am_cbps_1 <- bal.tab(intervention_1~confounders_1,
            weights=cbps_w,
            method="weighting",
            s.d.denom="pooled")$Balance$Diff.Adj

lpcbps1 <- love.plot(bal.tab(intervention_1~confounders_1_plot,
                             weights=cbps_w,
                             method="weighting",
                             s.d.denom="pooled"),
                     colors = c("grey", "black"),
                     var.order = "unadjusted",
                     shapes=c("circle","square"),
                     stars = "raw", var.names = v, legend="none", position=1, abs=TRUE)

surveyDesigncbps_1  <- svydesign(ids=1:n,weights=~cbps_w,data=data_imputed_1)
# s1<-svykm(Surv(TIME_CHD, failure)~totp, design=surveyDesign)
modelcbps_1   <- svycoxph(Surv(TIME_CHD, failure) ~ totp,design=surveyDesigncbps_1) 
summary(modelcbps_1)



timeCBPS2    <- tryCatch(  system.time( cbps_w  <- CBPS(as.factor(intervention_2) ~ scale(f45multi)+scale(f45mvmin)+scale(ethnic)+scale(parity)+scale(booph)+scale(meno)+scale(brca_f2)+scale(colon_f2)+scale(endo_f2)+scale(skin_f2)+scale(melan_f2)+scale(
  othca10y)+scale(dvt)+scale(stroke)+scale(mi)+scale(diab)+scale(hicholrp)+scale(osteopor)+scale(cvd)+
    scale(cabg)+scale(atrialfb)+scale(aortican)+scale(angina)+scale(hip55)+scale(smokevr)+scale(alcohol)+
    scale(fruits)+scale(vegtabls)+scale(f60enrgy)+scale(syst_bl)+scale(dias_bl)+scale(bmix_bl)+scale(educ)+
    scale(income),method = "exact", data=data_imputed_2, ATT=0)$weights )
  , error=function(e) NULL)#endtryCatch


n <- dim(confounders_2)[1]
am_cbps_2 <- bal.tab(intervention_2~confounders_2,
            weights=cbps_w,
            method="weighting",
            s.d.denom="pooled")$Balance$Diff.Adj

lpcbps2 <- love.plot(bal.tab(intervention_2~confounders_2_plot,
                             weights=cbps_w,
                             method="weighting",
                             s.d.denom="pooled"),
                     colors = c("grey", "black"),
                     var.order = "unadjusted",
                     shapes=c("circle","square"),
                     stars = "raw", var.names = v, legend="none", position=1, abs=TRUE)

surveyDesigncbps_2  <- svydesign(ids=1:n,weights=~cbps_w,data=data_imputed_2)
modelcbps_2   <- svycoxph(Surv(TIME_CHD, failure) ~ totp,design=surveyDesigncbps_2) 
summary(modelcbps_2)





timeCBPS3    <- tryCatch(  system.time( cbps_w  <- CBPS(as.factor(intervention_3) ~ scale(f45multi)+scale(f45mvmin)+scale(ethnic)+scale(parity)+scale(booph)+scale(meno)+scale(brca_f2)+scale(colon_f2)+scale(endo_f2)+scale(skin_f2)+scale(melan_f2)+scale(
  othca10y)+scale(dvt)+scale(stroke)+scale(mi)+scale(diab)+scale(hicholrp)+scale(osteopor)+scale(cvd)+
    scale(cabg)+scale(atrialfb)+scale(aortican)+scale(angina)+scale(hip55)+scale(smokevr)+scale(alcohol)+
    scale(fruits)+scale(vegtabls)+scale(f60enrgy)+scale(syst_bl)+scale(dias_bl)+scale(bmix_bl)+scale(educ)+
    scale(income),method = "exact", data=data_imputed_3, ATT=0)$weights )
  , error=function(e) NULL)#endtryCatch


n <- dim(confounders_3)[1]
am_cbps_3 <- bal.tab(intervention_3~confounders_3,
            weights=cbps_w,
            method="weighting",
            s.d.denom="pooled")$Balance$Diff.Adj

lpcbps3 <- love.plot(bal.tab(intervention_3~confounders_3_plot,
                             weights=cbps_w,
                             method="weighting",
                             s.d.denom="pooled"),
                     colors = c("grey", "black"),
                     var.order = "unadjusted",
                     shapes=c("circle","square"),
                     stars = "raw", var.names = v, legend="none", position=1, abs=TRUE)

surveyDesigncbps_3  <- svydesign(ids=1:n,weights=~cbps_w,data=data_imputed_3)
modelcbps_3   <- svycoxph(Surv(TIME_CHD, failure) ~ totp,design=surveyDesigncbps_3) 
summary(modelcbps_3)








####################################################################################################
####################################################################################################
####################################################################################################

#EBAL
timeEBAL1    <- tryCatch(  system.time(ebal1   <- weightit(totp~f45multi+f45mvmin+ethnic+parity+
                                                             booph+meno+brca_f2+colon_f2+endo_f2+skin_f2+melan_f2+
                                                             othca10y+dvt+stroke+mi+diab+hicholrp+osteopor+cvd+cabg+
                                                             atrialfb+aortican+angina+hip55+smokevr+alcohol+fruits+
                                                             vegtabls+f60enrgy+syst_bl+dias_bl+bmix_bl+educ+income , 
                                                           data = data_imputed_1,
                                                           method = "ebal", 
                                                           estimand = "ATE"))
                           , error=function(e) NULL)#endtryCatch


n <- dim(confounders_1)[1]

lpebal1 <- love.plot(bal.tab(intervention_1~confounders_1_plot,
                             weights=ebal1$weights,
                             method="weighting",
                             s.d.denom="pooled"),
                     colors = c("grey", "black"),
                     var.order = "unadjusted",
                     shapes=c("circle","square"),
                     stars = "raw")

am_eba_1 <- bal.tab(intervention_1~confounders_1,
                    weights=ebal1$weights,
                    method="weighting",
                    s.d.denom="pooled")$Balance$Diff.Adj

surveyDesignebal_1  <- svydesign(ids=1:n,weights=~ebal1$weights,data=data_imputed_1)
# s1<-svykm(Surv(TIME_CHD, failure)~totp, design=surveyDesign)
modelebal_1  <- svycoxph(Surv(TIME_CHD, failure) ~ totp,design=surveyDesignebal_1) 
summary(modelebal_1)



timeEBAL2    <- tryCatch(  system.time(ebal2   <- weightit(totp~f45multi+f45mvmin+ethnic+parity+
                                                             booph+meno+brca_f2+colon_f2+endo_f2+skin_f2+melan_f2+
                                                             othca10y+dvt+stroke+mi+diab+hicholrp+osteopor+cvd+cabg+
                                                             atrialfb+aortican+angina+hip55+smokevr+alcohol+fruits+
                                                             vegtabls+f60enrgy+syst_bl+dias_bl+bmix_bl+educ+income , 
                                                           data = data_imputed_2,
                                                           method = "ebal", 
                                                           estimand = "ATE"))
                           , error=function(e) NULL)#endtryCatch


n <- dim(confounders_2)[1]

lpebal2 <- love.plot(bal.tab(intervention_2~confounders_2_plot,
                             weights=ebal2$weights,
                             method="weighting",
                             s.d.denom="pooled"),
                     colors = c("grey", "black"),
                     var.order = "unadjusted",
                     shapes=c("circle","square"),
                     stars = "raw")

am_eba_2 <- bal.tab(intervention_2~confounders_2,
                    weights=ebal2$weights,
                    method="weighting",
                    s.d.denom="pooled")$Balance$Diff.Adj

surveyDesignebal_2  <- svydesign(ids=1:n,weights=~ebal2$weights,data=data_imputed_2)
modelebal_2  <- svycoxph(Surv(TIME_CHD, failure) ~ totp,design=surveyDesignebal_2) 
summary(modelebal_2)


timeEBAL3    <- tryCatch(  system.time(ebal3   <- weightit(totp~f45multi+f45mvmin+ethnic+parity+
                                                             booph+meno+brca_f2+colon_f2+endo_f2+skin_f2+melan_f2+
                                                             othca10y+dvt+stroke+mi+diab+hicholrp+osteopor+cvd+cabg+
                                                             atrialfb+aortican+angina+hip55+smokevr+alcohol+fruits+
                                                             vegtabls+f60enrgy+syst_bl+dias_bl+bmix_bl+educ+income , 
                                                           data = data_imputed_3,
                                                           method = "ebal", 
                                                           estimand = "ATE"))
                           , error=function(e) NULL)#endtryCatch


n <- dim(confounders_3)[1]

lpebal3 <- love.plot(bal.tab(intervention_3~confounders_3_plot,
                             weights=ebal3$weights,
                             method="weighting",
                             s.d.denom="pooled"),
                     colors = c("grey", "black"),
                     var.order = "unadjusted",
                     shapes=c("circle","square"),
                     stars = "raw")

am_eba_3 <- bal.tab(intervention_3~confounders_3,
                    weights=ebal3$weights,
                    method="weighting",
                    s.d.denom="pooled")$Balance$Diff.Adj

surveyDesignebal_3  <- svydesign(ids=1:n,weights=~ebal3$weights,data=data_imputed_3)
# s1<-svykm(Surv(TIME_CHD, failure)~totp, design=surveyDesign)
modelebal_3  <- svycoxph(Surv(TIME_CHD, failure) ~ totp,design=surveyDesignebal_3) 
summary(modelebal_3)


####################################################################################################
####################################################################################################
####################################################################################################
#PSM

sl.lib_ps  <- c('SL.glm','SL.glmnet',"SL.randomForest")

#Balancing covariates between treatment groups (binary)
Tr    <- abs(intervention_1-1)
timePSM1    <- tryCatch(  system.time( SL    <- SuperLearner(Tr, confounders_1 ,
                                                             SL.library = sl.lib_ps, 
                                                             family = "binomial"))
                          , error=function(e) NULL)#endtryCatch
ps    <- SL$SL.predict
ipw   <- ( Tr/ps + (1-Tr)/(1-ps)  )
summary(ipw)

timePSM1b    <- tryCatch(  system.time( m1 <- Match(Tr = Tr, X = ipw, estimand = "ATE")) , error=function(e) NULL)#endtryCatch

lpma1 <- love.plot(bal.tab(m1,intervention_1~confounders_1_plot),
                   colors = c("grey", "black"),
                   var.order = "unadjusted",
                   shapes=c("circle","square"),
                   stars = "raw",
                   var.names = v)
am_psm_1 <- bal.tab(m1,intervention_1~confounders_1)$Balance$Diff.Adj

matches <- factor(rep(m1$index.treated, 2))
matchedsample <- cbind(matches, data_imputed_1[c(m1$index.control, m1$index.treated),])
fitma1 <- coxph(Surv(TIME_CHD, failure) ~ totp + strata(matches), data=matchedsample)
summary(fitma1)



#Balancing covariates between treatment groups (binary)
Tr    <- abs(intervention_2-1)
timePSM2    <- tryCatch(  system.time( SL    <- SuperLearner(Tr, confounders_2 ,
                                                             SL.library = sl.lib_ps, 
                                                             family = "binomial"))
                          , error=function(e) NULL)#endtryCatch
ps    <- SL$SL.predict
ipw   <- ( Tr/ps + (1-Tr)/(1-ps)  )
summary(ipw)

timePSM2b    <- tryCatch(  system.time( m2 <- Match(Tr = Tr, X = ipw, estimand = "ATE")) , error=function(e) NULL)#endtryCatch


lpma2 <- love.plot(bal.tab(m2,intervention_2~confounders_2_plot),
                   colors = c("grey", "black"),
                   var.order = "unadjusted",
                   shapes=c("circle","square"),
                   stars = "raw",
                   var.names = v)
am_psm_2 <- bal.tab(m2,intervention_2~confounders_2)$Balance$Diff.Adj

matches <- factor(rep(m2$index.treated, 2))
matchedsample <- cbind(matches, data_imputed_2[c(m2$index.control, m2$index.treated),])
fitma2 <- coxph(Surv(TIME_CHD, failure) ~ totp + strata(matches), data=matchedsample)
summary(fitma2)


#Balancing covariates between treatment groups (binary)
Tr    <- abs(intervention_3-1)
timePSM3    <- tryCatch(  system.time( SL    <- SuperLearner(Tr, confounders_3 ,
                                                             SL.library = sl.lib_ps, 
                                                             family = "binomial"))
                          , error=function(e) NULL)#endtryCatch
ps    <- SL$SL.predict
ipw   <- ( Tr/ps + (1-Tr)/(1-ps)  )
summary(ipw)

timePSM3b    <- tryCatch(  system.time( m3 <- Match(Tr = Tr, X = ipw, estimand = "ATE")) , error=function(e) NULL)#endtryCatch


lpma3 <- love.plot(bal.tab(m3,intervention_3~confounders_3_plot),
                   colors = c("grey", "black"),
                   var.order = "unadjusted",
                   shapes=c("circle","square"),
                   stars = "raw",
                   var.names = v)
am_psm_3 <- bal.tab(m3,intervention_3~confounders_3)$Balance$Diff.Adj

matches <- factor(rep(m3$index.treated, 2))
matchedsample <- cbind(matches, data_imputed_3[c(m3$index.control, m3$index.treated),])
fitma3 <- coxph(Surv(TIME_CHD, failure) ~ totp + strata(matches), data=matchedsample)
summary(fitma3)


####################################################################################################
####################################################################################################
####################################################################################################
#GBM

stop.method <- "es.mean"

data_imputed_1$totp_01 <- abs(data_imputed_1$totp-1)

timegbm1    <- tryCatch( system.time(ps_gbm_1 <- weightit(totp_01~f45multi+f45mvmin+ethnic+parity+
                                                            booph+meno+brca_f2+colon_f2+endo_f2+skin_f2+melan_f2+
                                                            othca10y+dvt+stroke+mi+diab+hicholrp+osteopor+cvd+cabg+
                                                            atrialfb+aortican+angina+hip55+smokevr+alcohol+fruits+
                                                            vegtabls+f60enrgy+syst_bl+dias_bl+bmix_bl+educ+income,
                                                          data = data_imputed_1,
                                                          method = "gbm",
                                                          estimand = 'ATE',
                                                          stop.method = stop.method) )
                         , error=function(e) NULL)#endtryCatch

n <- dim(data_imputed_1)[1]
am_gbm_1 <- bal.tab(intervention_1~confounders_1,
            weights=ps_gbm_1$weights,
            method="weighting",
            s.d.denom="pooled")$Balance$Diff.Adj

lpgbm1 <- love.plot(bal.tab(intervention_1~confounders_1_plot,
                            weights=ps_gbm_1$weights,
                            method="weighting",
                            s.d.denom="pooled"),
                    colors = c("grey", "black"),
                    var.order = "unadjusted",
                    shapes=c("circle","square"),
                    stars = "raw")

surveyDesigngbm_1  <- svydesign(ids=1:n,weights=~ps_gbm_1$weights,data=data_imputed_1)
modelgbm_1   <- svycoxph(Surv(TIME_CHD, failure) ~ totp,design=surveyDesigngbm_1) 
summary(modelgbm_1)




data_imputed_2$totp_01 <- abs(data_imputed_2$totp-1)

timegbm2    <- tryCatch( system.time(ps_gbm_2 <- weightit(totp_01~f45multi+f45mvmin+ethnic+parity+
                                                            booph+meno+brca_f2+colon_f2+endo_f2+skin_f2+melan_f2+
                                                            othca10y+dvt+stroke+mi+diab+hicholrp+osteopor+cvd+cabg+
                                                            atrialfb+aortican+angina+hip55+smokevr+alcohol+fruits+
                                                            vegtabls+f60enrgy+syst_bl+dias_bl+bmix_bl+educ+income,
                                                          data = data_imputed_2,
                                                          method = "gbm",
                                                          estimand = 'ATE',
                                                          stop.method = stop.method) )
                         , error=function(e) NULL)#endtryCatch

n <- dim(data_imputed_2)[1]
am_gbm_2 <- bal.tab(intervention_2~confounders_2,
            weights=ps_gbm_2$weights,
            method="weighting",
            s.d.denom="pooled")$Balance$Diff.Adj

lpgbm2 <- love.plot(bal.tab(intervention_2~confounders_2_plot,
                            weights=ps_gbm_2$weights,
                            method="weighting",
                            s.d.denom="pooled"),
                    colors = c("grey", "black"),
                    var.order = "unadjusted",
                    shapes=c("circle","square"),
                    stars = "raw")

surveyDesigngbm_2  <- svydesign(ids=1:n,weights=~ps_gbm_2$weights,data=data_imputed_2)
modelgbm_2   <- svycoxph(Surv(TIME_CHD, failure) ~ totp,design=surveyDesigngbm_2) 
summary(modelgbm_2)




data_imputed_3$totp_01 <- abs(data_imputed_3$totp-1)

timegbm3    <- tryCatch( system.time(ps_gbm_3 <- weightit(totp_01~f45multi+f45mvmin+ethnic+parity+
                                                            booph+meno+brca_f2+colon_f2+endo_f2+skin_f2+melan_f2+
                                                            othca10y+dvt+stroke+mi+diab+hicholrp+osteopor+cvd+cabg+
                                                            atrialfb+aortican+angina+hip55+smokevr+alcohol+fruits+
                                                            vegtabls+f60enrgy+syst_bl+dias_bl+bmix_bl+educ+income,
                                                          data = data_imputed_3,
                                                          method = "gbm",
                                                          estimand = 'ATE',
                                                          stop.method = stop.method) )
                         , error=function(e) NULL)#endtryCatch

n <- dim(data_imputed_3)[1]
am_gbm_3 <- bal.tab(intervention_3~confounders_3,
            weights=ps_gbm_3$weights,
            method="weighting",
            s.d.denom="pooled")$Balance$Diff.Adj

lpgbm2 <- love.plot(bal.tab(intervention_3~confounders_3_plot,
                            weights=ps_gbm_3$weights,
                            method="weighting",
                            s.d.denom="pooled"),
                    colors = c("grey", "black"),
                    var.order = "unadjusted",
                    shapes=c("circle","square"),
                    stars = "raw")

surveyDesigngbm_3  <- svydesign(ids=1:n,weights=~ps_gbm_3$weights,data=data_imputed_3)
modelgbm_3   <- svycoxph(Surv(TIME_CHD, failure) ~ totp,design=surveyDesigngbm_3) 
summary(modelgbm_3)




####################################################################################################
####################################################################################################
####################################################################################################
#BalSL

sl.lib_ps  <- c('SL.glm','SL.glmnet',"SL.randomForest")

stop.method <- "es.mean"

data_imputed_1$totp_01 <- abs(data_imputed_1$totp-1)

timebalSL1    <- tryCatch( system.time(ps_balSL_1 <- weightit(totp_01~f45multi+f45mvmin+ethnic+parity+
                                                                booph+meno+brca_f2+colon_f2+endo_f2+skin_f2+melan_f2+
                                                                othca10y+dvt+stroke+mi+diab+hicholrp+osteopor+cvd+cabg+
                                                                atrialfb+aortican+angina+hip55+smokevr+alcohol+fruits+
                                                                vegtabls+f60enrgy+syst_bl+dias_bl+bmix_bl+educ+income,
                                                              data = data_imputed_1,
                                                              method = "super",
                                                              estimand = 'ATE',
                                                              SL.library =  sl.lib_ps ,
                                                              SL.method = "method.balance",
                                                              stop.method = stop.method) )
                           , error=function(e) NULL)#endtryCatch

n <- dim(data_imputed_1)[1]
am_bSL_1 <- bal.tab(intervention_1~confounders_1,
            weights=ps_balSL_1$weights,
            method="weighting",
            s.d.denom="pooled")$Balance$Diff.Adj

lpbalSL1 <- love.plot(bal.tab(intervention_1~confounders_1_plot,
                              weights=ps_balSL_1$weights,
                              method="weighting",
                              s.d.denom="pooled"),
                      colors = c("grey", "black"),
                      var.order = "unadjusted",
                      shapes=c("circle","square"),
                      stars = "raw")

surveyDesignbalSL_1  <- svydesign(ids=1:n,weights=~ps_balSL_1$weights,data=data_imputed_1)
modelbalSL_1   <- svycoxph(Surv(TIME_CHD, failure) ~ totp,design=surveyDesignbalSL_1) 
summary(modelbalSL_1)




data_imputed_2$totp_01 <- abs(data_imputed_2$totp-1)

timebalSL2    <- tryCatch( system.time(ps_balSL_2 <- weightit(totp_01~f45multi+f45mvmin+ethnic+parity+
                                                                booph+meno+brca_f2+colon_f2+endo_f2+skin_f2+melan_f2+
                                                                othca10y+dvt+stroke+mi+diab+hicholrp+osteopor+cvd+cabg+
                                                                atrialfb+aortican+angina+hip55+smokevr+alcohol+fruits+
                                                                vegtabls+f60enrgy+syst_bl+dias_bl+bmix_bl+educ+income,
                                                              data = data_imputed_2,
                                                              method = "super",
                                                              estimand = 'ATE',
                                                              SL.library =  sl.lib_ps ,
                                                              SL.method = "method.balance",
                                                              stop.method = stop.method) )
                           , error=function(e) NULL)#endtryCatch

n <- dim(data_imputed_2)[1]
am_bSL_2 <- bal.tab(intervention_2~confounders_2,
            weights=ps_balSL_2$weights,
            method="weighting",
            s.d.denom="pooled")$Balance$Diff.Adj

lpbalSL2 <- love.plot(bal.tab(intervention_2~confounders_2_plot,
                              weights=ps_balSL_2$weights,
                              method="weighting",
                              s.d.denom="pooled"),
                      colors = c("grey", "black"),
                      var.order = "unadjusted",
                      shapes=c("circle","square"),
                      stars = "raw")

surveyDesignbalSL_2  <- svydesign(ids=1:n,weights=~ps_balSL_2$weights,data=data_imputed_2)
modelbalSL_2   <- svycoxph(Surv(TIME_CHD, failure) ~ totp,design=surveyDesignbalSL_2) 
summary(modelbalSL_2)




data_imputed_3$totp_01 <- abs(data_imputed_3$totp-1)

timebalSL3    <- tryCatch( system.time(ps_balSL_3 <- weightit(totp_01~f45multi+f45mvmin+ethnic+parity+
                                                                booph+meno+brca_f2+colon_f2+endo_f2+skin_f2+melan_f2+
                                                                othca10y+dvt+stroke+mi+diab+hicholrp+osteopor+cvd+cabg+
                                                                atrialfb+aortican+angina+hip55+smokevr+alcohol+fruits+
                                                                vegtabls+f60enrgy+syst_bl+dias_bl+bmix_bl+educ+income,
                                                              data = data_imputed_3,
                                                              method = "super",
                                                              estimand = 'ATE',
                                                              SL.library =  sl.lib_ps ,
                                                              SL.method = "method.balance",
                                                              stop.method = stop.method) )
                           , error=function(e) NULL)#endtryCatch

n <- dim(data_imputed_3)[1]
am_bSL_3 <- bal.tab(intervention_3~confounders_3,
            weights=ps_balSL_3$weights,
            method="weighting",
            s.d.denom="pooled")$Balance$Diff.Adj

lpbalSL3 <- love.plot(bal.tab(intervention_3~confounders_3_plot,
                              weights=ps_balSL_3$weights,
                              method="weighting",
                              s.d.denom="pooled"),
                      colors = c("grey", "black"),
                      var.order = "unadjusted",
                      shapes=c("circle","square"),
                      stars = "raw")

surveyDesignbalSL_3  <- svydesign(ids=1:n,weights=~ps_balSL_3$weights,data=data_imputed_3)
modelbalSL_3   <- svycoxph(Surv(TIME_CHD, failure) ~ totp,design=surveyDesignbalSL_3) 
summary(modelbalSL_3)

####################################################################################################
####################################################################################################
####################################################################################################

#SBW
data_sbw_1      <- data_imputed_1
data_sbw_1$totp      <- abs(data_sbw_1$totp-1)
balsbw         <-  list()
balsbw$bal_cov <- c("f45multi","f45mvmin","ethnic","parity","booph","meno","brca_f2","colon_f2","endo_f2","skin_f2","melan_f2",
                    "othca10y","dvt","stroke","mi","diab","hicholrp","osteopor","cvd",
                    "cabg","atrialfb","aortican","angina","hip55","smokevr","alcohol",
                    "fruits","vegtabls","f60enrgy","syst_bl","dias_bl","bmix_bl","educ","income")
balsbw$bal_tol <- 0.0001
balsbw$bal_std <- TRUE
balsbw$bal_alg <- FALSE
balsbw$bal_gri <- 0.1
balsbw$bal_sam <- 1

timesbw1    <- tryCatch( system.time(out1 <-  sbw(data_sbw_1, ind = "totp",
             sol = list(sol_nam = "gurobi"),
             bal = balsbw,
             par = list(par_est = "ate")
             ) )
             , error=function(e) NULL)#endtryCatch

sbw <- out1$dat_weights$weights

n <- dim(confounders_1)[1]
am_sbw_1 <- bal.tab(intervention_1~confounders_1,
            weights=sbw,
            method="weighting",
            s.d.denom="pooled")$Balance$Diff.Adj

lpsbw1 <- love.plot(bal.tab(intervention_1~confounders_1_plot,
                  weights=sbw,
                  method="weighting",
                  s.d.denom="pooled"),
          colors = c("grey", "black"),
          var.order = "unadjusted",
          shapes=c("circle","square"),
          stars = "raw")

surveyDesignsbw_1  <- svydesign(ids=1:n,weights=~sbw,data=data_imputed_1)
modelsbw_1   <- svycoxph(Surv(TIME_CHD, failure) ~ totp,design=surveyDesignsbw_1)
summary(modelsbw_1)



###
data_sbw_2     <- data_imputed_2
data_sbw_2$totp      <- abs(data_sbw_2$totp-1)
timesbw2    <- tryCatch( system.time(out2 <-  sbw(data_sbw_2, ind = "totp",
                                                  sol = list(sol_nam = "gurobi"),
                                                  bal = balsbw,
                                                  par = list(par_est = "ate")
) )
, error=function(e) NULL)#endtryCatch

sbw <- out2$dat_weights$weights

n <- dim(confounders_2)[1]
am_sbw_2 <- bal.tab(intervention_2~confounders_2,
            weights=sbw,
            method="weighting",
            s.d.denom="pooled")$Balance$Diff.Adj

lpsbw2 <- love.plot(bal.tab(intervention_2~confounders_2_plot,
                  weights=sbw,
                  method="weighting",
                  s.d.denom="pooled"),
          colors = c("grey", "black"),
          var.order = "unadjusted",
          shapes=c("circle","square"),
          stars = "raw")

surveyDesignsbw_2  <- svydesign(ids=1:n,weights=~sbw,data=data_imputed_2)
modelsbw_2   <- svycoxph(Surv(TIME_CHD, failure) ~ totp,design=surveyDesignsbw_2)
summary(modelsbw_2)


###
data_sbw_3     <- data_imputed_3
data_sbw_3$totp      <- abs(data_sbw_3$totp-1)
timesbw3    <- tryCatch( system.time(out3 <-  sbw(data_sbw_3, ind = "totp",
             sol = list(sol_nam = "gurobi"),
             bal = balsbw,
             par = list(par_est = "ate")
))
, error=function(e) NULL)#endtryCatch
#
sbw <- out3$dat_weights$weights

n <- dim(confounders_3)[1]
am_sbw_3 <- bal.tab(intervention_3~confounders_3,
            weights=sbw,
            method="weighting",
            s.d.denom="pooled")$Balance$Diff.Adj

lpsbw3 <- love.plot(bal.tab(intervention_3~confounders_3_plot,
                  weights=sbw,
                  method="weighting",
                  s.d.denom="pooled"),
          colors = c("grey", "black"),
          var.order = "unadjusted",
          shapes=c("circle","square"),
          stars = "raw")

surveyDesignsbw_3  <- svydesign(ids=1:n,weights=~sbw,data=data_imputed_3)
modelsbw_3   <- svycoxph(Surv(TIME_CHD, failure) ~ totp,design=surveyDesignsbw_3)
summary(modelsbw_3)

am_una_1 <- bal.tab(intervention_1~confounders_1)$Balance$Diff.Un
am_una_2 <- bal.tab(intervention_2~confounders_2)$Balance$Diff.Un
am_una_3 <- bal.tab(intervention_3~confounders_3)$Balance$Diff.Un

save.image("case_study_b2.Rdata")    

####################################################################################################
####################################################################################################
####################################################################################################
#Tables

load("case_study_b2.Rdata")    

table_1 <- rbind(c(exp(modelw_1$coefficients[1]),exp(confint(modelw_1)), sqrt(vcov(modelw_1)),
                   timeROW1[1]+timeROW1[2],median(abs(am_row_1)), min(abs(am_row_1)), max(abs(am_row_1)) ),
                 c(exp(modelipw_1$coefficients[1]),exp(confint(modelipw_1)), sqrt(vcov(modelipw_1)),
                   timeIPW1[1]+timeIPW1[2] ,median(abs(am_ipw_1)), min(abs(am_ipw_1)), max(abs(am_ipw_1))),
                 c(exp(modelbalSL_1$coefficients[1]),exp(confint(modelbalSL_1)), sqrt(vcov(modelbalSL_1)),
                   timebalSL1[1]+timebalSL1[2] ,median(abs(am_bSL_1)), min(abs(am_bSL_1)), max(abs(am_bSL_1))),
                 c(exp(modelgbm_1$coefficients[1]),exp(confint(modelgbm_1)), sqrt(vcov(modelgbm_1)),
                   timegbm1[1]+timegbm1[2] ,median(abs(am_gbm_1)), min(abs(am_gbm_1)), max(abs(am_gbm_1))),
                 c(exp(modelcbps_1$coefficients[1]),exp(confint(modelcbps_1)), sqrt(vcov(modelcbps_1)),
                   timeCBPS1[1]+timeCBPS1[2] ,median(abs(am_cbps_1)), min(abs(am_cbps_1)), max(abs(am_cbps_1))),
                 c(exp(modelebal_1$coefficients[1]),exp(confint(modelebal_1)), sqrt(vcov(modelebal_1)),
                   timeEBAL1[1]+timeEBAL1[2] ,median(abs(am_eba_1)), min(abs(am_eba_1)), max(abs(am_eba_1))),
                 c(exp(fitma1$coefficients[1]),exp(confint(fitma1)), sqrt(vcov(fitma1)),
                   timePSM1[1]+timePSM1[2] + timePSM1b[1]+timePSM1b[2] ,median(abs(am_psm_1)), min(abs(am_psm_1)), max(abs(am_psm_1))),
                 c(exp(modelm_1$coefficients[1]),exp(confint(modelm_1))[1,] , sqrt(vcov(modelm_1)[1,1]),
                   timeOM1[1]+timeOM1[2] ,c(0,0,0)),
                 c(exp(modeln_1$coefficients[1]),exp(confint(modeln_1)), sqrt(vcov(modeln_1)),
                   timeNa1[1]+timeNa1[2] ,median(abs(am_una_1)), min(abs(am_una_1)), max(abs(am_una_1)))
                 )

rownames(table_1) <- c("ROW","IPW","BalSL","GBM","CBPS","EBAL","PSM","OM","Naive")
colnames(table_1) <- c("HR","2.5%","97.5%","SE","Time in seconds", "Median absolute balance", "Min absolute balance", "Max absolute balance")

table_1

table_2 <- rbind(c(exp(modelw_2$coefficients[1]),exp(confint(modelw_2)), sqrt(vcov(modelw_2)),
                   timeROW2[1]+timeROW2[2],median(abs(am_row_2)), min(abs(am_row_2)), max(abs(am_row_2)) ),
                 c(exp(modelipw_2$coefficients[1]),exp(confint(modelipw_2)), sqrt(vcov(modelipw_2)),
                   timeIPW2[1]+timeIPW2[2] ,median(abs(am_ipw_2)), min(abs(am_ipw_2)), max(abs(am_ipw_2))),
                 c(exp(modelbalSL_2$coefficients[1]),exp(confint(modelbalSL_2)), sqrt(vcov(modelbalSL_2)),
                   timebalSL2[1]+timebalSL2[2] ,median(abs(am_bSL_2)), min(abs(am_bSL_2)), max(abs(am_bSL_2))),
                 c(exp(modelgbm_2$coefficients[1]),exp(confint(modelgbm_2)), sqrt(vcov(modelgbm_2)),
                   timegbm2[1]+timegbm2[2] ,median(abs(am_gbm_2)), min(abs(am_gbm_2)), max(abs(am_gbm_2))),
                 c(exp(modelcbps_2$coefficients[1]),exp(confint(modelcbps_2)), sqrt(vcov(modelcbps_2)),
                   timeCBPS2[1]+timeCBPS2[2] ,median(abs(am_cbps_2)), min(abs(am_cbps_2)), max(abs(am_cbps_2))),
                 c(exp(modelebal_2$coefficients[1]),exp(confint(modelebal_2)), sqrt(vcov(modelebal_2)),
                   timeEBAL2[1]+timeEBAL2[2] ,median(abs(am_eba_2)), min(abs(am_eba_2)), max(abs(am_eba_2))),
                 c(exp(fitma2$coefficients[1]),exp(confint(fitma2)), sqrt(vcov(fitma2)),
                   timePSM2[1]+timePSM2[2] + timePSM2b[1]+timePSM2b[2] ,median(abs(am_psm_2)), min(abs(am_psm_2)), max(abs(am_psm_2))),
                 c(exp(modelm_2$coefficients[1]),exp(confint(modelm_2))[1,] , sqrt(vcov(modelm_2)[1,1]),
                   timeOM2[1]+timeOM2[2] ,c(0,0,0)),
                 c(exp(modeln_2$coefficients[1]),exp(confint(modeln_2)), sqrt(vcov(modeln_2)),
                   timeNa2[1]+timeNa2[2] ,median(abs(am_una_2)), min(abs(am_una_2)), max(abs(am_una_2)))
)

rownames(table_2) <- c("ROW","IPW","BalSL","GBM","CBPS","EBAL","PSM","OM","Naive")
colnames(table_2) <- c("HR","2.5%","97.5%","SE","Time in seconds", "Median absolute balance", "Min absolute balance", "Max absolute balance")

table_2

table_3 <- rbind(c(exp(modelw_3$coefficients[1]),exp(confint(modelw_3)), sqrt(vcov(modelw_3)),
                   timeROW3[1]+timeROW3[2],median(abs(am_row_3)), min(abs(am_row_3)), max(abs(am_row_3)) ),
                 c(exp(modelipw_3$coefficients[1]),exp(confint(modelipw_3)), sqrt(vcov(modelipw_3)),
                   timeIPW3[1]+timeIPW3[2] ,median(abs(am_ipw_3)), min(abs(am_ipw_3)), max(abs(am_ipw_3))),
                 c(exp(modelbalSL_3$coefficients[1]),exp(confint(modelbalSL_3)), sqrt(vcov(modelbalSL_3)),
                   timebalSL3[1]+timebalSL3[2] ,median(abs(am_bSL_3)), min(abs(am_bSL_3)), max(abs(am_bSL_3))),
                 c(exp(modelgbm_3$coefficients[1]),exp(confint(modelgbm_3)), sqrt(vcov(modelgbm_3)),
                   timegbm3[1]+timegbm3[2] ,median(abs(am_gbm_3)), min(abs(am_gbm_3)), max(abs(am_gbm_3))),
                 c(exp(modelcbps_3$coefficients[1]),exp(confint(modelcbps_3)), sqrt(vcov(modelcbps_3)),
                   timeCBPS3[1]+timeCBPS3[2] ,median(abs(am_cbps_3)), min(abs(am_cbps_3)), max(abs(am_cbps_3))),
                 c(exp(modelebal_3$coefficients[1]),exp(confint(modelebal_3)), sqrt(vcov(modelebal_3)),
                   timeEBAL3[1]+timeEBAL3[2] ,median(abs(am_eba_3)), min(abs(am_eba_3)), max(abs(am_eba_3))),
                 c(exp(fitma2$coefficients[1]),exp(confint(fitma2)), sqrt(vcov(fitma2)),
                   timePSM3[1]+timePSM3[2] + timePSM3b[1]+timePSM3b[2] ,median(abs(am_psm_3)), min(abs(am_psm_3)), max(abs(am_psm_3))),
                 c(exp(modelm_3$coefficients[1]),exp(confint(modelm_3))[1,] , sqrt(vcov(modelm_3)[1,1]),
                   timeOM3[1]+timeOM3[2] ,c(0,0,0)),
                 c(exp(modeln_3$coefficients[1]),exp(confint(modeln_3)), sqrt(vcov(modeln_3)),
                   timeNa3[1]+timeNa3[2] ,median(abs(am_una_3)), min(abs(am_una_3)), max(abs(am_una_3)))
)

rownames(table_3) <- c("ROW","IPW","BalSL","GBM","CBPS","EBAL","PSM","OM","Naive")
colnames(table_3) <- c("HR","2.5%","97.5%","Time in seconds", "Median absolute balance", "Min absolute balance", "Max absolute balance")

table_3

write.csv(table_1,"table_1.csv")
write.csv(table_2,"table_2.csv")
write.csv(table_3,"table_3.csv")





#########################################################################################################
grid.newpage()
pushViewport(viewport(layout = grid.layout(2, 2, heights = unit(c(0.5, 5, 
                                                                  0.5,  5), "null"))))  

grid.text("Time since menopause: 0-10 Years", 
          vp = viewport(layout.pos.row = 1, layout.pos.col = 1:2),gp=gpar(fontsize=20))


print(lprow1, vp = viewport(layout.pos.row = 2, layout.pos.col = 1))
print(km1$plot, vp = viewport(layout.pos.row = 2, layout.pos.col = 2))


#########################################################################################################
grid.newpage()
pushViewport(viewport(layout = grid.layout(2, 2, heights = unit(c(0.5, 5, 
                                                                  0.5,  5), "null"))))  

grid.text("Time since menopause: 11-20 Years", 
          vp = viewport(layout.pos.row = 1, layout.pos.col = 1:2),gp=gpar(fontsize=20))


print(lprow2, vp = viewport(layout.pos.row = 2, layout.pos.col = 1))
print(km2$plot, vp = viewport(layout.pos.row = 2, layout.pos.col = 2))


#########################################################################################################
grid.newpage()
pushViewport(viewport(layout = grid.layout(2, 2, heights = unit(c(0.5, 5, 
                                                                  0.5,  5), "null"))))  

grid.text("Time since menopause: 20+ Years", 
          vp = viewport(layout.pos.row = 1, layout.pos.col = 1:2),gp=gpar(fontsize=20))


print(lprow3, vp = viewport(layout.pos.row = 2, layout.pos.col = 1))
print(km3$plot, vp = viewport(layout.pos.row = 2, layout.pos.col = 2))


#########################################################################################################
#######KM
grid.newpage()
pushViewport(viewport(layout = grid.layout(4, 2, heights = unit(c(0.5, 5, 0.5,  5), "null"))))  

grid.text("Time since menopause: 0-10 Years", 
          vp = viewport(layout.pos.row = 1, layout.pos.col = 1),gp=gpar(fontsize=14))

grid.text("Time since menopause: 11-20 Years", 
          vp = viewport(layout.pos.row = 1, layout.pos.col = 2),gp=gpar(fontsize=14))

print(km1$plot, vp = viewport(layout.pos.row = 2, layout.pos.col = 1))
print(km2$plot, vp = viewport(layout.pos.row = 2, layout.pos.col = 2))

grid.text("Time since menopause: 20+ Years", 
          vp = viewport(layout.pos.row = 3, layout.pos.col = 1),gp=gpar(fontsize=14))

print(km3$plot, vp = viewport(layout.pos.row = 4, layout.pos.col = 1))


