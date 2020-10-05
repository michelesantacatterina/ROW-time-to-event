########################################################################
########################################################################
#***********************************************************************
#
#   In this file, we provide code to make all the plots in the original
#   manuscript and in the supplementary material. 
# 
#***********************************************************************
########################################################################
########################################################################


library(ggplot2)
library(ggpubr)
library(grid)

rm(list = ls())



n <- 1000

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


      #########################################################################################################
      #########################################################################################################
      # Positivity
      #########################################################################################################
      #########################################################################################################
      
      SIZEL   <- 2
      SIZET   <- 14
      SIZELe  <- 14
      LENGTHLe <- 2.2
      
      #Notes: Binary: 
      # * While lack of overlap increases, 
      #   - IPW keeps low bias but performs worse in terms of MSE
      # * While misspecification increases:
      #   - bias OM and IPW increases
      
      load(paste("data_results/simu_n",n,"positivity.Rdata"))
      
      #########################################################################################################
      #######BIAS
      dfBIAS_b <- dfBIAS[which(dfBIAS$Treatment=='binary'),]
      dfBIAS_c <- dfBIAS[which(dfBIAS$Treatment=='continuous'),]
      
      #grey  <- c("grey80", "grey75", "grey70", "grey60", "grey50", "grey30", "grey20", "grey0")
      
      #           IPW       BalSL     GBM         PSM         CBPS      EBAL      SBW         OM          naive     ROW  
      grey  <- c("#0066CC", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#990000", "#000000")
      #           IPW      BalSL     GBM         PSM         CBPS      EBAL      SBW         OM      naive     ROW   
      lines <- c("dashed","dashed","dashed",   "dotted",   "dotdash","dotdash","dotdash", "dotted", "dotted", "solid")
      
      
      gBiasP_b <- ggplot(dfBIAS_b, aes(x=Positivity,y=abs(Bias), col=Method)) + 
        geom_line(aes(linetype=Method),size=SIZEL) +
        coord_cartesian(ylim = c(0, 0.5)) +
        xlab("Strenght of PPV") + 
        scale_x_continuous(breaks=c(0.1, 1, 2),labels=c("Weak", "Moderate", "Strong")) +
        theme_bw() + theme(legend.position="none", 
                           axis.title = element_text(size = SIZET),
                           axis.text = element_text(size = SIZET),
                           axis.text.x = element_text(hjust = c(0.5,0.5,0.8)) ) +
        ylab("Abs(Bias)") + 
        scale_color_manual(values=grey) +
        scale_linetype_manual(values = lines) +
        scale_fill_manual(values=grey,
                          name="Method",
                          breaks=unique(dfBIAS_b$Method),
                          labels=unique(dfBIAS_b$Method))
      
      #             IPW     BalSL       GBM         CBPS          OM        naive       ROW       npCBPS
      grey  <- c("#0066CC", "#E69F00", "#56B4E9",  "#F0E442",  "#CC79A7",  "#990000", "#000000", "#0072B2")
      #             IPW     BalSL       GBM         CBPS          OM        naive       ROW    npCBPS
      lines <- c("dashed","dashed",  "dashed",   "dotdash",    "dotted",   "dotted",  "solid", "dotdash") 
      
      
      gBiasP_c <- ggplot(dfBIAS_c, aes(x=Positivity,y=abs(Bias), col=Method)) + 
        geom_line(aes(linetype=Method),size=SIZEL) +
        coord_cartesian(ylim = c(0, 5)) +
        xlab("Strenght of PPV") + 
        scale_x_continuous(breaks=c(1, 3, 5),labels=c("Weak", "Moderate", "Strong")) + 
        theme_bw() + theme(legend.position="none", 
                           axis.title = element_text(size = SIZET),
                           axis.text = element_text(size = SIZET),
                           axis.text.x = element_text(hjust = c(0.5,0.5,0.8)) ) +
        ylab("Abs(Bias)") + 
        scale_color_manual(values=grey) +
        scale_linetype_manual(values = lines) +
        scale_fill_manual(values=grey,
                          name="Method",
                          breaks=unique(dfBIAS_c$Method),
                          labels=unique(dfBIAS_c$Method))
      
      
      #########################################################################################################
      #######MSE
      dfMSE_b <- dfMSE[which(dfMSE$Treatment=='binary'),]
      dfMSE_c <- dfMSE[which(dfMSE$Treatment=='continuous'),]
      
      
      #           IPW       BalSL     GBM         PSM         CBPS      EBAL      SBW         OM          naive     ROW  
      grey  <- c("#0066CC", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#990000", "#000000")
      #           IPW      BalSL     GBM         PSM         CBPS      EBAL      SBW         OM      naive     ROW   
      lines <- c("dashed","dashed","dashed",   "dotted",   "dotdash","dotdash","dotdash", "dotted", "dotted", "solid")
      
      
      gMseP_b <- ggplot(dfMSE_b, aes(x=Positivity,y=sqrt(MSE), col=Method)) + 
        geom_line(aes(linetype=Method),size=SIZEL) +
        coord_cartesian(ylim = c(0, 0.5)) +
        xlab("Strenght of PPV") + 
        scale_x_continuous(breaks=c(0.1, 1, 2),labels=c("Weak", "Moderate", "Strong")) +
        theme_bw() + theme(legend.position="none", 
                           axis.title = element_text(size = SIZET),
                           axis.text = element_text(size = SIZET),
                           axis.text.x = element_text(hjust = c(0.5,0.5,0.8)) ) +
        ylab("RMSE")+ 
        scale_color_manual(values=grey) +
        scale_linetype_manual(values = lines) +
        scale_fill_manual(values=grey,
                          name="Method",
                          breaks=unique(dfMSE_b$Method),
                          labels=unique(dfMSE_b$Method))
      
      
      #             IPW     BalSL       GBM         CBPS          OM        naive       ROW       npCBPS
      grey  <- c("#0066CC", "#E69F00", "#56B4E9",  "#F0E442",  "#CC79A7",  "#990000", "#000000", "#0072B2")
      #             IPW     BalSL       GBM         CBPS          OM        naive       ROW    npCBPS
      lines <- c("dashed","dashed",  "dashed",   "dotdash",    "dotted",   "dotted",  "solid", "dotdash") 
      
      
      gMseP_c <- ggplot(dfMSE_c, aes(x=Positivity,y=sqrt(MSE), col=Method)) + 
        geom_line(aes(linetype=Method),size=SIZEL) +
        coord_cartesian(ylim = c(0, 5)) +
        xlab("Strenght of PPV") + 
        scale_x_continuous(breaks=c(1, 3, 5),labels=c("Weak", "Moderate", "Strong")) + 
        theme_bw() + theme(legend.position="none", 
                           axis.title = element_text(size = SIZET),
                           axis.text = element_text(size = SIZET),
                           axis.text.x = element_text(hjust = c(0.5,0.5,0.8)) ) +
        ylab("RMSE")+ 
        scale_color_manual(values=grey) +
        scale_linetype_manual(values = lines) +
        scale_fill_manual(values=grey,
                          name="Method",
                          breaks=unique(dfMSE_c$Method),
                          labels=unique(dfMSE_c$Method))
      
      
      #########################################################################################################
      #######BALANCE
      dfBAL_b <- dfBAL[which(dfBAL$Treatment=='binary'),]
      dfBAL_c <- dfBAL[which(dfBAL$Treatment=='continuous'),]
      
      
      #           IPW       BalSL     GBM         PSM         CBPS      EBAL      SBW         OM          naive     ROW  
      grey  <- c("#0066CC", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#990000", "#000000")
      #           IPW      BalSL     GBM         PSM         CBPS      EBAL      SBW         OM      naive     ROW   
      lines <- c("dashed","dashed","dashed",   "dotted",   "dotdash","dotdash","dotdash", "dotted", "dotted", "solid")
      
      
      gBalP_b <- ggplot(dfBAL_b, aes(x=Positivity,y=Balance, col=Method)) + 
        geom_line(aes(linetype=Method),size=SIZEL) +
        coord_cartesian(ylim = c(0, 1)) +
        xlab("Strenght of PPV") + 
        scale_x_continuous(breaks=c(0.1, 1, 2),labels=c("Weak", "Moderate", "Strong")) +
        theme_bw() + theme(legend.position="none", 
                           axis.title = element_text(size = SIZET),
                           axis.text = element_text(size = SIZET),
                           axis.text.x = element_text(hjust = c(0.5,0.5,0.8)) ) +
        ylab("Balance")+ 
        scale_color_manual(values=grey) +
        scale_linetype_manual(values = lines) +
        scale_fill_manual(values=grey,
                          name="Method",
                          breaks=unique(dfBAL_b$Method),
                          labels=unique(dfBAL_b$Method))
      
      
      #             IPW     BalSL       GBM         CBPS          OM        naive       ROW       npCBPS
      grey  <- c("#0066CC", "#E69F00", "#56B4E9",  "#F0E442",  "#CC79A7",  "#990000", "#000000", "#0072B2")
      #             IPW     BalSL       GBM         CBPS          OM        naive       ROW    npCBPS
      lines <- c("dashed","dashed",  "dashed",   "dotdash",    "dotted",   "dotted",  "solid", "dotdash") 
      
      
      
      gBalP_c <- ggplot(dfBAL_c, aes(x=Positivity,y=Balance, col=Method)) + 
        geom_line(aes(linetype=Method),size=SIZEL) +
        coord_cartesian(ylim = c(0, 1)) +
        xlab("Strenght of PPV") + 
        scale_x_continuous(breaks=c(1, 3, 5),labels=c("Weak", "Moderate", "Strong")) + 
        theme_bw() + theme(legend.position="none", 
                           axis.title = element_text(size = SIZET),
                           axis.text = element_text(size = SIZET),
                           axis.text.x = element_text(hjust = c(0.5,0.5,0.8)) ) +
        ylab("Balance")+ 
        scale_color_manual(values=grey) +
        scale_linetype_manual(values = lines) +
        scale_fill_manual(values=grey,
                          name="Method",
                          breaks=unique(dfBAL_c$Method),
                          labels=unique(dfBAL_c$Method))
      
      
      
      #########################################################################################################
      #######TIME
      dfTIME_b <- dfTIME[which(dfTIME$Treatment=='binary'),]
      dfTIME_c <- dfTIME[which(dfTIME$Treatment=='continuous'),]
      
      
      #           IPW       BalSL     GBM         PSM         CBPS      EBAL      SBW         OM          naive     ROW  
      grey  <- c("#0066CC", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#990000", "#000000")
      #           IPW      BalSL     GBM         PSM         CBPS      EBAL      SBW         OM      naive     ROW   
      lines <- c("dashed","dashed","dashed",   "dotted",   "dotdash","dotdash","dotdash", "dotted", "dotted", "solid")
      
      
      gTimeP_b <- ggplot(dfTIME_b, aes(x=Positivity,y=Time, col=Method)) + 
        geom_line(aes(linetype=Method),size=SIZEL) +
        xlab("Strenght of PPV") + 
        scale_x_continuous(breaks=c(0.1, 1, 2),labels=c("Weak", "Moderate", "Strong")) + 
        theme_bw() + theme(legend.position="none", 
                           axis.title = element_text(size = SIZET),
                           axis.text = element_text(size = SIZET),
                           axis.text.x = element_text(hjust = c(0.5,0.5,0.8)) ) +
        ylab("Time")+ 
        scale_color_manual(values=grey) +
        scale_linetype_manual(values = lines) +
        scale_fill_manual(values=grey,
                          name="Method",
                          breaks=unique(dfTIME_b$Method),
                          labels=unique(dfTIME_b$Method))
      
      
      #             IPW     BalSL       GBM         CBPS          OM        naive       ROW       npCBPS
      grey  <- c("#0066CC", "#E69F00", "#56B4E9",  "#F0E442",  "#CC79A7",  "#990000", "#000000", "#0072B2")
      #             IPW     BalSL       GBM         CBPS          OM        naive       ROW    npCBPS
      lines <- c("dashed","dashed",  "dashed",   "dotdash",    "dotted",   "dotted",  "solid", "dotdash") 
      
      
      
      gTimeP_c <- ggplot(dfTIME_c, aes(x=Positivity,y=Time, col=Method)) + 
        geom_line(aes(linetype=Method),size=SIZEL) +
        xlab("Strenght of PPV") + 
        scale_x_continuous(breaks=c(1, 3, 5),labels=c("Weak", "Moderate", "Strong")) + 
        theme_bw() + theme(legend.position="none", 
                           axis.title = element_text(size = SIZET),
                           axis.text = element_text(size = SIZET),
                           axis.text.x = element_text(hjust = c(0.5,0.5,0.8)) ) +
        ylab("Time")+ 
        scale_color_manual(values=grey) +
        scale_linetype_manual(values = lines) +
        scale_fill_manual(values=grey,
                          name="Method",
                          breaks=unique(dfTIME_c$Method),
                          labels=unique(dfTIME_c$Method))
      
      
      #########################################################################################################
      #########################################################################################################
      # Misspecification
      #########################################################################################################
      #########################################################################################################
      
      
      load(paste("data_results/simu_n",n,"misspecified.Rdata"))
      
      
      #########################################################################################################
      #######BIAS
      dfBIAS_b <- dfBIAS[which(dfBIAS$Treatment=='binary'),]
      dfBIAS_c <- dfBIAS[which(dfBIAS$Treatment=='continuous'),]
      
      
      
      #           IPW       BalSL     GBM         PSM         CBPS      EBAL      SBW         OM          naive     ROW  
      grey  <- c("#0066CC", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#990000", "#000000")
      #           IPW      BalSL     GBM         PSM         CBPS      EBAL      SBW         OM      naive     ROW   
      lines <- c("dashed","dashed","dashed",   "dotted",   "dotdash","dotdash","dotdash", "dotted", "dotted", "solid")
      
      
      gBiasM_b <- ggplot(dfBIAS_b, aes(x=Misspecification,y=abs(Bias), col=Method)) + 
        geom_line(aes(linetype=Method),size=SIZEL) +
        coord_cartesian(ylim = c(0, 1)) +
        xlab("Strenght of misspecification") + 
        scale_x_continuous(breaks=c(0, 0.5, 1),labels=c("Null", "Moderate", "Strong")) + 
        theme_bw() + theme(legend.position="none", 
                           axis.title = element_text(size = SIZET),
                           axis.text = element_text(size = SIZET),
                           axis.text.x = element_text(hjust = c(0.5,0.5,0.8)) ) +
        ylab("Abs(Bias)") + 
        scale_color_manual(values=grey) +
        scale_linetype_manual(values = lines) +
        scale_fill_manual(values=grey,
                          name="Method",
                          breaks=unique(dfBIAS_b$Method),
                          labels=unique(dfBIAS_b$Method))
      
      
      #             IPW     BalSL       GBM         CBPS          OM        naive       ROW       npCBPS
      grey  <- c("#0066CC", "#E69F00", "#56B4E9",  "#F0E442",  "#CC79A7",  "#990000", "#000000", "#0072B2")
      #             IPW     BalSL       GBM         CBPS          OM        naive       ROW    npCBPS
      lines <- c("dashed","dashed",  "dashed",   "dotdash",    "dotted",   "dotted",  "solid", "dotdash") 
      
      
      gBiasM_c <- ggplot(dfBIAS_c, aes(x=Misspecification,y=abs(Bias), col=Method)) + 
        geom_line(aes(linetype=Method),size=SIZEL) +
        coord_cartesian(ylim = c(0, 0.3)) +
        xlab("Strenght of misspecification") + 
        scale_x_continuous(breaks=c(0, 0.5, 1),labels=c("Null", "Moderate", "Strong")) + 
        theme_bw() + theme(legend.position="none", 
                           axis.title = element_text(size = SIZET),
                           axis.text = element_text(size = SIZET),
                           axis.text.x = element_text(hjust = c(0.5,0.5,0.8)) ) +
        ylab("Abs(Bias)") + 
        scale_color_manual(values=grey) +
        scale_linetype_manual(values = lines) +
        scale_fill_manual(values=grey,
                          name="Method",
                          breaks=unique(dfBIAS_c$Method),
                          labels=unique(dfBIAS_c$Method))
      
      
      #########################################################################################################
      #######MSE
      dfMSE_b <- dfMSE[which(dfMSE$Treatment=='binary'),]
      dfMSE_c <- dfMSE[which(dfMSE$Treatment=='continuous'),]
      
      
      
      #           IPW       BalSL     GBM         PSM         CBPS      EBAL      SBW         OM          naive     ROW  
      grey  <- c("#0066CC", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#990000", "#000000")
      #           IPW      BalSL     GBM         PSM         CBPS      EBAL      SBW         OM      naive     ROW   
      lines <- c("dashed","dashed","dashed",   "dotted",   "dotdash","dotdash","dotdash", "dotted", "dotted", "solid")
      
      
      gMseM_b <- ggplot(dfMSE_b, aes(x=Misspecification,y=sqrt(MSE), col=Method)) + 
        geom_line(aes(linetype=Method),size=SIZEL) +
        xlab("Strenght of misspecification") + 
        scale_x_continuous(breaks=c(0, 0.5, 1),labels=c("Null", "Moderate", "Strong")) + 
        coord_cartesian(ylim = c(0, 1)) +
        theme_bw() + theme(legend.position="none", 
                           axis.title = element_text(size = SIZET),
                           axis.text = element_text(size = SIZET),
                           axis.text.x = element_text(hjust = c(0.5,0.5,0.8)) ) +
        ylab("RMSE")+ 
        scale_color_manual(values=grey) +
        scale_linetype_manual(values = lines) +
        scale_fill_manual(values=grey,
                          name="Method",
                          breaks=unique(dfMSE_b$Method),
                          labels=unique(dfMSE_b$Method))
      
      
      #             IPW     BalSL       GBM         CBPS          OM        naive       ROW       npCBPS
      grey  <- c("#0066CC", "#E69F00", "#56B4E9",  "#F0E442",  "#CC79A7",  "#990000", "#000000", "#0072B2")
      #             IPW     BalSL       GBM         CBPS          OM        naive       ROW    npCBPS
      lines <- c("dashed","dashed",  "dashed",   "dotdash",    "dotted",   "dotted",  "solid", "dotdash") 
      
      
      
      gMseM_c <- ggplot(dfMSE_c, aes(x=Misspecification,y=sqrt(MSE), col=Method)) + 
        geom_line(aes(linetype=Method),size=SIZEL) +
        coord_cartesian(ylim = c(0, 0.3)) +
        xlab("Strenght of misspecification") + 
        scale_x_continuous(breaks=c(0, 0.5, 1),labels=c("Null", "Moderate", "Strong")) + 
        theme_bw() + theme(legend.position="none", 
                           axis.title = element_text(size = SIZET),
                           axis.text = element_text(size = SIZET),
                           axis.text.x = element_text(hjust = c(0.5,0.5,0.8)) ) +
        ylab("RMSE")+ 
        scale_color_manual(values=grey) +
        scale_linetype_manual(values = lines) +
        scale_fill_manual(values=grey,
                          name="Method",
                          breaks=unique(dfMSE_c$Method),
                          labels=unique(dfMSE_c$Method))
      
      
      #########################################################################################################
      #######BALANCE
      dfBAL_b <- dfBAL[which(dfBAL$Treatment=='binary'),]
      dfBAL_c <- dfBAL[which(dfBAL$Treatment=='continuous'),]
      
      
      #           IPW       BalSL     GBM         PSM         CBPS      EBAL      SBW         OM          naive     ROW  
      grey  <- c("#0066CC", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#990000", "#000000")
      #           IPW      BalSL     GBM         PSM         CBPS      EBAL      SBW         OM      naive     ROW   
      lines <- c("dashed","dashed","dashed",   "dotted",   "dotdash","dotdash","dotdash", "dotted", "dotted", "solid")
      
      
      gBalM_b <- ggplot(dfBAL_b, aes(x=Misspecification,y=Balance, col=Method)) + 
        geom_line(aes(linetype=Method),size=SIZEL) +
        xlab("Strenght of misspecification") + 
        scale_x_continuous(breaks=c(0, 0.5, 1),labels=c("Null", "Moderate", "Strong")) + 
        theme_bw() + theme(legend.position="none", 
                           axis.title = element_text(size = SIZET),
                           axis.text = element_text(size = SIZET),
                           axis.text.x = element_text(hjust = c(0.5,0.5,0.8)) ) +
        ylab("Balance")+ 
        scale_color_manual(values=grey) +
        scale_linetype_manual(values = lines) +
        scale_fill_manual(values=grey,
                          name="Method",
                          breaks=unique(dfBAL_b$Method),
                          labels=unique(dfBAL_b$Method))
      
      
      #             IPW     BalSL       GBM         CBPS          OM        naive       ROW       npCBPS
      grey  <- c("#0066CC", "#E69F00", "#56B4E9",  "#F0E442",  "#CC79A7",  "#990000", "#000000", "#0072B2")
      #             IPW     BalSL       GBM         CBPS          OM        naive       ROW    npCBPS
      lines <- c("dashed","dashed",  "dashed",   "dotdash",    "dotted",   "dotted",  "solid", "dotdash") 
      
      
      gBalM_c <- ggplot(dfBAL_c, aes(x=Misspecification,y=Balance, col=Method)) + 
        geom_line(aes(linetype=Method),size=SIZEL) +
        xlab("Strenght of misspecification") + 
        scale_x_continuous(breaks=c(0, 0.5, 1),labels=c("Null", "Moderate", "Strong")) +  
        theme_bw() + theme(legend.position="none", 
                           axis.title = element_text(size = SIZET),
                           axis.text = element_text(size = SIZET),
                           axis.text.x = element_text(hjust = c(0.5,0.5,0.8)) ) +
        ylab("Balance")+ 
        scale_color_manual(values=grey) +
        scale_linetype_manual(values = lines) +
        scale_fill_manual(values=grey,
                          name="Method",
                          breaks=unique(dfBAL_c$Method),
                          labels=unique(dfBAL_c$Method))
      
      
      
      
      #########################################################################################################
      #######TIME
      dfTIME_b <- dfTIME[which(dfTIME$Treatment=='binary'),]
      dfTIME_c <- dfTIME[which(dfTIME$Treatment=='continuous'),]
      
      
      #           IPW       BalSL     GBM         PSM         CBPS      EBAL      SBW         OM          naive     ROW  
      grey  <- c("#0066CC", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#990000", "#000000")
      #           IPW      BalSL     GBM         PSM         CBPS      EBAL      SBW         OM      naive     ROW   
      lines <- c("dashed","dashed","dashed",   "dotted",   "dotdash","dotdash","dotdash", "dotted", "dotted", "solid")
      
      
      gTimeM_b <- ggplot(dfTIME_b, aes(x=Misspecification,y=Time, col=Method)) + 
        geom_line(aes(linetype=Method),size=SIZEL) +
        xlab("Strenght of misspecification") + 
        scale_x_continuous(breaks=c(0, 0.5, 1),labels=c("Null", "Moderate", "Strong")) + 
        theme_bw() + theme(legend.position="none", 
                           axis.title = element_text(size = SIZET),
                           axis.text = element_text(size = SIZET),
                           axis.text.x = element_text(hjust = c(0.5,0.5,0.8)) ) +
        ylab("Time")+ 
        scale_color_manual(values=grey) +
        scale_linetype_manual(values = lines) +
        scale_fill_manual(values=grey,
                          name="Method",
                          breaks=unique(dfTIME_b$Method),
                          labels=unique(dfTIME_b$Method))
      
      
      #             IPW     BalSL       GBM         CBPS          OM        naive       ROW       npCBPS
      grey  <- c("#0066CC", "#E69F00", "#56B4E9",  "#F0E442",  "#CC79A7",  "#990000", "#000000", "#0072B2")
      #             IPW     BalSL       GBM         CBPS          OM        naive       ROW    npCBPS
      lines <- c("dashed","dashed",  "dashed",   "dotdash",    "dotted",   "dotted",  "solid", "dotdash") 
      
      
      
      gTimeM_c <- ggplot(dfTIME_c, aes(x=Misspecification,y=Time, col=Method)) + 
        geom_line(aes(linetype=Method),size=SIZEL) +
        xlab("Strenght of misspecification") + 
        scale_x_continuous(breaks=c(0, 0.5, 1),labels=c("Null", "Moderate", "Strong")) + 
        theme_bw() + theme(legend.position="none", 
                           axis.title = element_text(size = SIZET),
                           axis.text = element_text(size = SIZET),
                           axis.text.x = element_text(hjust = c(0.5,0.5,0.8)) ) +
        ylab("Time")+ 
        scale_color_manual(values=grey) +
        scale_linetype_manual(values = lines) +
        scale_fill_manual(values=grey,
                          name="Method",
                          breaks=unique(dfTIME_c$Method),
                          labels=unique(dfTIME_c$Method))
      
      
      
      
      
      
      
      
      
      #########################################################################################################
      #########################################################################################################
      # Censoring
      #########################################################################################################
      #########################################################################################################
      
      
      load(paste("data_results/simu_n",n,"censoring.Rdata"))
      
      
      #########################################################################################################
      #######BIAS
      dfBIAS_b <- dfBIAS[which(dfBIAS$Treatment=='binary'),]
      dfBIAS_c <- dfBIAS[which(dfBIAS$Treatment=='continuous'),]
      
      
      #           IPW       BalSL     GBM         PSM         CBPS      EBAL      SBW         OM          naive     ROW  
      grey  <- c("#0066CC", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#990000", "#000000")
      #           IPW      BalSL     GBM         PSM         CBPS      EBAL      SBW         OM      naive     ROW   
      lines <- c("dashed","dashed","dashed",   "dotted",   "dotdash","dotdash","dotdash", "dotted", "dotted", "solid")
      
      
      gBiasC_b <- ggplot(dfBIAS_b, aes(x=Censoring,y=abs(Bias), col=Method)) + 
        geom_line(aes(linetype=Method),size=SIZEL) +
        coord_cartesian(ylim = c(0, 1)) +
        xlab("Percentage of censoring") + 
        scale_x_continuous(breaks=c(0, 0.4, 0.8),labels=c("0%", "40%", "80%")) + 
        theme_bw() + theme(legend.position="none", 
                           axis.title = element_text(size = SIZET),
                           axis.text = element_text(size = SIZET),
                           axis.text.x = element_text(hjust = c(0.5,0.5,0.8)) ) +
        ylab("Abs(Bias)") + 
        scale_color_manual(values=grey) +
        scale_linetype_manual(values = lines) +
        scale_fill_manual(values=grey,
                          name="Method",
                          breaks=unique(dfBIAS_b$Method),
                          labels=unique(dfBIAS_b$Method))
      
      
      #             IPW     BalSL       GBM         CBPS          OM        naive       ROW       npCBPS
      grey  <- c("#0066CC", "#E69F00", "#56B4E9",  "#F0E442",  "#CC79A7",  "#990000", "#000000", "#0072B2")
      #             IPW     BalSL       GBM         CBPS          OM        naive       ROW    npCBPS
      lines <- c("dashed","dashed",  "dashed",   "dotdash",    "dotted",   "dotted",  "solid", "dotdash") 
      
      
      gBiasC_c <- ggplot(dfBIAS_c, aes(x=Censoring,y=abs(Bias), col=Method)) + 
        geom_line(aes(linetype=Method),size=SIZEL) +
        coord_cartesian(ylim = c(0, 0.5)) +
        xlab("Percentage of censoring") + 
        scale_x_continuous(breaks=c(0, 0.4, 0.8),labels=c("0%", "40%", "80%")) +
        theme_bw() + theme(legend.position="none", 
                           axis.title = element_text(size = SIZET),
                           axis.text = element_text(size = SIZET),
                           axis.text.x = element_text(hjust = c(0.5,0.5,0.8)) ) +
        ylab("Abs(Bias)") + 
        scale_color_manual(values=grey) +
        scale_linetype_manual(values = lines) +
        scale_fill_manual(values=grey,
                          name="Method",
                          breaks=unique(dfBIAS_c$Method),
                          labels=unique(dfBIAS_c$Method))
      
      
      #########################################################################################################
      #######MSE
      dfMSE_b <- dfMSE[which(dfMSE$Treatment=='binary'),]
      dfMSE_c <- dfMSE[which(dfMSE$Treatment=='continuous'),]
      
      
      
      #           IPW       BalSL     GBM         PSM         CBPS      EBAL      SBW         OM          naive     ROW  
      grey  <- c("#0066CC", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#990000", "#000000")
      #           IPW      BalSL     GBM         PSM         CBPS      EBAL      SBW         OM      naive     ROW   
      lines <- c("dashed","dashed","dashed",   "dotted",   "dotdash","dotdash","dotdash", "dotted", "dotted", "solid")
      
      
      gMseC_b <- ggplot(dfMSE_b, aes(x=Censoring,y=sqrt(MSE), col=Method)) + 
        geom_line(aes(linetype=Method),size=SIZEL) +
        xlab("Percentage of censoring") + 
        scale_x_continuous(breaks=c(0, 0.4, 0.8),labels=c("0%", "40%", "80%")) +
        coord_cartesian(ylim = c(0, 1)) +
        theme_bw() + theme(legend.position="none", 
                           axis.title = element_text(size = SIZET),
                           axis.text = element_text(size = SIZET),
                           axis.text.x = element_text(hjust = c(0.5,0.5,0.8)) ) +
        ylab("RMSE")+ 
        scale_color_manual(values=grey) +
        scale_linetype_manual(values = lines) +
        scale_fill_manual(values=grey,
                          name="Method",
                          breaks=unique(dfMSE_b$Method),
                          labels=unique(dfMSE_b$Method))
      
      
      #             IPW     BalSL       GBM         CBPS          OM        naive       ROW       npCBPS
      grey  <- c("#0066CC", "#E69F00", "#56B4E9",  "#F0E442",  "#CC79A7",  "#990000", "#000000", "#0072B2")
      #             IPW     BalSL       GBM         CBPS          OM        naive       ROW    npCBPS
      lines <- c("dashed","dashed",  "dashed",   "dotdash",    "dotted",   "dotted",  "solid", "dotdash") 
      
      
      
      gMseC_c <- ggplot(dfMSE_c, aes(x=Censoring,y=sqrt(MSE), col=Method)) + 
        geom_line(aes(linetype=Method),size=SIZEL) +
        coord_cartesian(ylim = c(0, 0.5)) +
        xlab("Percentage of censoring") + 
        scale_x_continuous(breaks=c(0, 0.4, 0.8),labels=c("0%", "40%", "80%")) +
        theme_bw() + theme(legend.position="none", 
                           axis.title = element_text(size = SIZET),
                           axis.text = element_text(size = SIZET),
                           axis.text.x = element_text(hjust = c(0.5,0.5,0.8)) ) +
        ylab("RMSE")+ 
        scale_color_manual(values=grey) +
        scale_linetype_manual(values = lines) +
        scale_fill_manual(values=grey,
                          name="Method",
                          breaks=unique(dfMSE_c$Method),
                          labels=unique(dfMSE_c$Method))
      
      
      #########################################################################################################
      #######BALANCE
      dfBAL_b <- dfBAL[which(dfBAL$Treatment=='binary'),]
      dfBAL_c <- dfBAL[which(dfBAL$Treatment=='continuous'),]
      
      
      #           IPW       BalSL     GBM         PSM         CBPS      EBAL      SBW         OM          naive     ROW  
      grey  <- c("#0066CC", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#990000", "#000000")
      #           IPW      BalSL     GBM         PSM         CBPS      EBAL      SBW         OM      naive     ROW   
      lines <- c("dashed","dashed","dashed",   "dotted",   "dotdash","dotdash","dotdash", "dotted", "dotted", "solid")
      
      
      gBalC_b <- ggplot(dfBAL_b, aes(x=Censoring,y=Balance, col=Method)) + 
        geom_line(aes(linetype=Method),size=SIZEL) +
        xlab("Percentage of censoring") + 
        scale_x_continuous(breaks=c(0, 0.4, 0.8),labels=c("0%", "40%", "80%")) +
        theme_bw() + theme(legend.position="none", 
                           axis.title = element_text(size = SIZET),
                           axis.text = element_text(size = SIZET),
                           axis.text.x = element_text(hjust = c(0.5,0.5,0.8)) ) +
        ylab("Balance")+ 
        scale_color_manual(values=grey) +
        scale_linetype_manual(values = lines) +
        scale_fill_manual(values=grey,
                          name="Method",
                          breaks=unique(dfBAL_b$Method),
                          labels=unique(dfBAL_b$Method))
      
      
      #             IPW     BalSL       GBM         CBPS          OM        naive       ROW       npCBPS
      grey  <- c("#0066CC", "#E69F00", "#56B4E9",  "#F0E442",  "#CC79A7",  "#990000", "#000000", "#0072B2")
      #             IPW     BalSL       GBM         CBPS          OM        naive       ROW    npCBPS
      lines <- c("dashed","dashed",  "dashed",   "dotdash",    "dotted",   "dotted",  "solid", "dotdash") 
      
      
      gBalC_c <- ggplot(dfBAL_c, aes(x=Censoring,y=Balance, col=Method)) + 
        geom_line(aes(linetype=Method),size=SIZEL) +
        xlab("Percentage of censoring") + 
        scale_x_continuous(breaks=c(0, 0.4, 0.8),labels=c("0%", "40%", "80%")) +
        theme_bw() + theme(legend.position="none", 
                           axis.title = element_text(size = SIZET),
                           axis.text = element_text(size = SIZET),
                           axis.text.x = element_text(hjust = c(0.5,0.5,0.8)) ) +
        ylab("Balance")+ 
        scale_color_manual(values=grey) +
        scale_linetype_manual(values = lines) +
        scale_fill_manual(values=grey,
                          name="Method",
                          breaks=unique(dfBAL_c$Method),
                          labels=unique(dfBAL_c$Method))
      
      
      
      
      #########################################################################################################
      #######TIME
      dfTIME_b <- dfTIME[which(dfTIME$Treatment=='binary'),]
      dfTIME_c <- dfTIME[which(dfTIME$Treatment=='continuous'),]
      
      
      #           IPW       BalSL     GBM         PSM         CBPS      EBAL      SBW         OM          naive     ROW  
      grey  <- c("#0066CC", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#990000", "#000000")
      #           IPW      BalSL     GBM         PSM         CBPS      EBAL      SBW         OM      naive     ROW   
      lines <- c("dashed","dashed","dashed",   "dotted",   "dotdash","dotdash","dotdash", "dotted", "dotted", "solid")
      
      
      gTimeC_b <- ggplot(dfTIME_b, aes(x=Censoring,y=Time, col=Method)) + 
        geom_line(aes(linetype=Method),size=SIZEL) +
        xlab("Percentage of censoring") + 
        scale_x_continuous(breaks=c(0, 0.4, 0.8),labels=c("0%", "40%", "80%")) +
        theme_bw() + theme(legend.position="none", 
                           axis.title = element_text(size = SIZET),
                           axis.text = element_text(size = SIZET),
                           axis.text.x = element_text(hjust = c(0.5,0.5,0.8)) ) +
        ylab("Time")+ 
        scale_color_manual(values=grey) +
        scale_linetype_manual(values = lines) +
        scale_fill_manual(values=grey,
                          name="Method",
                          breaks=unique(dfTIME_b$Method),
                          labels=unique(dfTIME_b$Method))
      
      
      
      
      Noprint_b <- ggplot(dfBAL_b, aes(x=Censoring,y=Balance, col=Method)) + 
        geom_line(aes(linetype=Method),size=SIZEL) +
        xlab("Percentage of censoring") + 
        scale_x_continuous(breaks=c(0, 0.4, 0.8),labels=c("0%", "40%", "80%")) +
        theme_bw() + theme(legend.position="bottom", 
                           axis.title = element_text(size = SIZET),
                           axis.text = element_text(size = SIZET),
                           legend.text = element_text(size = SIZELe),
                           legend.title = element_text(size = SIZELe),
                           legend.key.width = unit(LENGTHLe,"cm"),
                           axis.text.x = element_text(hjust = c(0.5,0.5,0.8)) ) +
        ylab("Balance")+ 
        scale_color_manual(values=grey) +
        scale_linetype_manual(values = lines) +
        scale_fill_manual(values=grey,
                          name="Method",
                          breaks=unique(dfBAL_b$Method),
                          labels=unique(dfBAL_b$Method))
      
      
      #             IPW     BalSL       GBM         CBPS          OM        naive       ROW       npCBPS
      grey  <- c("#0066CC", "#E69F00", "#56B4E9",  "#F0E442",  "#CC79A7",  "#990000", "#000000", "#0072B2")
      #             IPW     BalSL       GBM         CBPS          OM        naive       ROW    npCBPS
      lines <- c("dashed","dashed",  "dashed",   "dotdash",    "dotted",   "dotted",  "solid", "dotdash") 
      
      
      
      gTimeC_c <- ggplot(dfTIME_c, aes(x=Censoring,y=Time, col=Method)) + 
        geom_line(aes(linetype=Method),size=SIZEL) +
        xlab("Percentage of censoring") + 
        scale_x_continuous(breaks=c(0, 0.4, 0.8),labels=c("0%", "40%", "80%")) +
        theme_bw() + theme(legend.position="none", 
                           axis.title = element_text(size = SIZET),
                           axis.text = element_text(size = SIZET),
                           axis.text.x = element_text(hjust = c(0.5,0.5,0.8)) ) +
        ylab("Time")+ 
        scale_color_manual(values=grey) +
        scale_linetype_manual(values = lines) +
        scale_fill_manual(values=grey,
                          name="Method",
                          breaks=unique(dfTIME_c$Method),
                          labels=unique(dfTIME_c$Method))
      
      
      
      
    
      
      
      Noprint_c <- ggplot(dfTIME_c, aes(x=Censoring,y=Time, col=Method)) + 
        geom_line(aes(linetype=Method),size=SIZEL) +
        xlab("Percentage of censoring") + 
        scale_x_continuous(breaks=c(0, 0.4, 0.8),labels=c("0%", "40%", "80%")) +
        theme_bw() + theme(legend.position="bottom", 
                           axis.title = element_text(size = SIZET),
                           axis.text = element_text(size = SIZET),
                           legend.text = element_text(size = SIZELe),
                           legend.title = element_text(size = SIZELe),
                           legend.key.width = unit(LENGTHLe,"cm"),
                           axis.text.x = element_text(hjust = c(0.5,0.5,0.8)) ) +
        ylab("Time")+ 
        scale_color_manual(values=grey) +
        scale_linetype_manual(values = lines) +
        scale_fill_manual(values=grey,
                          name="Method",
                          breaks=unique(dfTIME_c$Method),
                          labels=unique(dfTIME_c$Method))
      
      
      
      
      #########################################################################################################
      #######BINARY TRT - Bias, MSE across positivity, misspecification and censoring
      # Extract the legend. Returns a gtable
      leg_b <- get_legend(Noprint_b)
      
      
      grid.newpage()
      pushViewport(viewport(layout = grid.layout(10, 2, heights = unit(c(0.5, 0.5, 5, 
                                                                        0.5, 0.5, 5, 
                                                                        0.5, 0.5, 5, 
                                                                        1), "null"))))  
      
      grid.text("Practical positivity violation (PPV)", 
                vp = viewport(layout.pos.row = 1, layout.pos.col = 1:2),gp=gpar(fontsize=20))
      
      
      print(gBiasP_b, vp = viewport(layout.pos.row = 3, layout.pos.col = 1))
      print(gMseP_b, vp = viewport(layout.pos.row = 3, layout.pos.col = 2))
      
      grid.text("Misspecification", vp = viewport(layout.pos.row = 4, layout.pos.col = 1:2),gp=gpar(fontsize=20))
      
      print(gBiasM_b, vp = viewport(layout.pos.row = 6, layout.pos.col = 1))
      print(gMseM_b, vp = viewport(layout.pos.row = 6, layout.pos.col = 2))
      
      grid.text("Censoring", vp = viewport(layout.pos.row = 7, layout.pos.col = 1:2),gp=gpar(fontsize=20))
      
      print(gBiasC_b, vp = viewport(layout.pos.row = 9, layout.pos.col = 1))
      print(gMseC_b, vp = viewport(layout.pos.row = 9, layout.pos.col = 2))
      
      print(as_ggplot(leg_b), vp = viewport(layout.pos.row = 10, layout.pos.col = 1:2),gp=gpar(fontsize=20))
      
      
      
      
      
      #########################################################################################################
      #######CONTIUNUOUS TRT - Bias, MSE across positivity, misspecification and censoring
      # Extract the legend. Returns a gtable
      leg_c <- get_legend(Noprint_c)
      
      
      grid.newpage()
      pushViewport(viewport(layout = grid.layout(10, 2, heights = unit(c(0.5, 0.5, 5, 
                                                                         0.5, 0.5, 5, 
                                                                         0.5, 0.5, 5, 
                                                                         1), "null"))))  
      
      grid.text("Practical positivity violation (PPV)", 
                vp = viewport(layout.pos.row = 1, layout.pos.col = 1:2),gp=gpar(fontsize=20))
      
      
      print(gBiasP_c, vp = viewport(layout.pos.row = 3, layout.pos.col = 1))
      print(gMseP_c, vp = viewport(layout.pos.row = 3, layout.pos.col = 2))
      
      grid.text("Misspecification", vp = viewport(layout.pos.row = 4, layout.pos.col = 1:2),gp=gpar(fontsize=20))
      
      print(gBiasM_c, vp = viewport(layout.pos.row = 6, layout.pos.col = 1))
      print(gMseM_c, vp = viewport(layout.pos.row = 6, layout.pos.col = 2))
      
      grid.text("Censoring", vp = viewport(layout.pos.row = 7, layout.pos.col = 1:2),gp=gpar(fontsize=20))
      
      print(gBiasC_c, vp = viewport(layout.pos.row = 9, layout.pos.col = 1))
      print(gMseC_c, vp = viewport(layout.pos.row = 9, layout.pos.col = 2))
      
      print(as_ggplot(leg_c), vp = viewport(layout.pos.row = 10, layout.pos.col = 1:2))
      
      
      
      
      
      
      
      
      #########################################################################################################
      #######BINARY TRT - Balance and Time across positivity, misspecification and censoring
      # Extract the legend. Returns a gtable
      leg_b <- get_legend(Noprint_b)
      
      
      grid.newpage()
      pushViewport(viewport(layout = grid.layout(10, 2, heights = unit(c(0.5, 0.5, 5, 
                                                                         0.5, 0.5, 5, 
                                                                         0.5, 0.5, 5, 
                                                                         1), "null"))))  
      
      grid.text("Practical positivity violation", 
                vp = viewport(layout.pos.row = 1, layout.pos.col = 1:2))
      
      
      print(gBalP_b, vp = viewport(layout.pos.row = 3, layout.pos.col = 1))
      print(gTimeP_b, vp = viewport(layout.pos.row = 3, layout.pos.col = 2))
      
      grid.text("Misspecification", vp = viewport(layout.pos.row = 4, layout.pos.col = 1:2))
      
      print(gBalM_b, vp = viewport(layout.pos.row = 6, layout.pos.col = 1))
      print(gTimeM_b, vp = viewport(layout.pos.row = 6, layout.pos.col = 2))
      
      grid.text("Censoring", vp = viewport(layout.pos.row = 7, layout.pos.col = 1:2))
      
      print(gBalC_b, vp = viewport(layout.pos.row = 9, layout.pos.col = 1))
      print(gTimeC_b, vp = viewport(layout.pos.row = 9, layout.pos.col = 2))
      
      print(as_ggplot(leg_b), vp = viewport(layout.pos.row = 10, layout.pos.col = 1:2))
      
      
      
      
      
      
      #########################################################################################################
      #######CONTINUOUS TRT - Balance and Time across positivity, misspecification and censoring
      # Extract the legend. Returns a gtable
      leg_c <- get_legend(Noprint_c)
      
      
      grid.newpage()
      pushViewport(viewport(layout = grid.layout(10, 2, heights = unit(c(0.5, 0.5, 5, 
                                                                         0.5, 0.5, 5, 
                                                                         0.5, 0.5, 5, 
                                                                         1), "null"))))  
      
      grid.text("Practical positivity violation", 
                vp = viewport(layout.pos.row = 1, layout.pos.col = 1:2))
      
      
      print(gBalP_c, vp = viewport(layout.pos.row = 3, layout.pos.col = 1))
      print(gTimeP_c, vp = viewport(layout.pos.row = 3, layout.pos.col = 2))
      
      grid.text("Misspecification", vp = viewport(layout.pos.row = 4, layout.pos.col = 1:2))
      
      print(gBalM_c, vp = viewport(layout.pos.row = 6, layout.pos.col = 1))
      print(gTimeM_c, vp = viewport(layout.pos.row = 6, layout.pos.col = 2))
      
      grid.text("Censoring", vp = viewport(layout.pos.row = 7, layout.pos.col = 1:2))
      
      print(gBalC_c, vp = viewport(layout.pos.row = 9, layout.pos.col = 1))
      print(gTimeC_c, vp = viewport(layout.pos.row = 9, layout.pos.col = 2))
      
      print(as_ggplot(leg_c), vp = viewport(layout.pos.row = 10, layout.pos.col = 1:2))
      
      
      # ggarrange(gBiasP_b,gMseP_b,gBiasM_b,gMseM_b,gBiasC_b,gMseC_b,nrow=3,ncol=2)
      # ggarrange(gBiasP_c,gMseP_c,gBiasM_c,gMseM_c,gBiasC_c,gMseC_c,nrow=3,ncol=2)
      # 
      # ggarrange(gBalP_b,gTimeP_b,gBalM_b,gTimeM_b,gBalC_b,gTimeC_b,nrow=3,ncol=2)
      # ggarrange(gBalP_c,gTimeP_c,gBalM_c,gTimeM_c,gBalC_c,gTimeC_c,nrow=3,ncol=2)
      
      # ggarrange(gBiasP_b,gBiasM_b,gBiasP_c,gBiasM_c,nrow=2,ncol=2)
      # ggarrange(gMseP_b,gMseM_b,gMseP_c,gMseM_c,nrow=2,ncol=2)
      # ggarrange(gBalP_b,gBalM_b,gBalP_c,gBalM_c,nrow=2,ncol=2)
      # ggarrange(gTimeP_b,gTimeM_b,gTimeP_c,gTimeM_c,nrow=2,ncol=2)


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

      n <- 1000
      
      SIZEL   <- 2
      SIZET   <- 14
      SIZELe  <- 14
      LENGTHLe <- 2.2
      
      #########################################################################################################
      #########################################################################################################
      # Positivity
      #########################################################################################################
      #########################################################################################################
      
      load(paste("data_results/simu_errors_n",n,"positivity.Rdata"))
      
      #########################################################################################################
      #######SE
      dfSE_b <- dfSE[which(dfSE$Treatment=='binary'),]
      dfSE_c <- dfSE[which(dfSE$Treatment=='continuous'),]
      
      #             Boot  Empirical   Naive   Robust
      grey  <- c("grey30", "grey0", "grey70", "grey50")
      #             Boot  Empirical   Naive   Robust
      lines <- c("dashed","solid","dotted","dotdash") 
      
      dfSE_b$Temp <- dfSE_b$SE[which(dfSE_b$Type=="Empirical")]/dfSE_b$SE
      gSEP_b <- ggplot(dfSE_b, aes(x=Positivity,y=Temp, col=Type)) + 
        geom_line(aes(linetype=Type),size=1) +
        coord_cartesian(ylim = c(0, 1)) +
        xlab("Strenght of practical positivity violation") + 
        scale_x_continuous(breaks=c(0.1, 1, 2),labels=c("Weak", "Moderate", "Strong")) + 
        theme_bw() + theme(legend.position="none", 
                           axis.title = element_text(size = SIZET),
                           axis.text = element_text(size = SIZET),
                           axis.text.x = element_text(hjust = c(0.5,0.5,0.8)) ) +
        ylab("SD/SE Ratio") + 
        scale_color_manual(values=c(grey[1],grey[2],grey[3],grey[4])) +
        scale_linetype_manual(values = lines) +
        scale_fill_manual(values=c(grey[1],grey[2],grey[3],grey[4]),
                          name="Method",
                          breaks=unique(dfSE_b$Method),
                          labels=unique(dfSE_b$Method))
      
      dfSE_c$Temp <- dfSE_c$SE[which(dfSE_c$Type=="Empirical")]/dfSE_c$SE
      gSEP_c <- ggplot(dfSE_c, aes(x=Positivity,y=Temp, col=Type)) + 
        geom_line(aes(linetype=Type),size=1) +
        coord_cartesian(ylim = c(0, 1.2)) +
        xlab("Strenght of practical positivity violation") + 
        scale_x_continuous(breaks=c(1, 3, 5),labels=c("Weak", "Moderate", "Strong")) + 
        theme_bw() + theme(legend.position="none", 
                           axis.title = element_text(size = SIZET),
                           axis.text = element_text(size = SIZET),
                           axis.text.x = element_text(hjust = c(0.5,0.5,0.8)) ) +
        ylab("SD/SE Ratio") + 
        scale_color_manual(values=c(grey[1],grey[2],grey[3],grey[4])) +
        scale_linetype_manual(values = lines) +
        scale_fill_manual(values=c(grey[1],grey[2],grey[3],grey[4]),
                          name="Method",
                          breaks=unique(dfSE_c$Method),
                          labels=unique(dfSE_c$Method))
      
      
      #########################################################################################################
      #######COVERAGE
      dfCO_b <- dfCO[which(dfCO$Treatment=='binary'),]
      dfCO_c <- dfCO[which(dfCO$Treatment=='continuous'),]
      
      
      #             Boot   Naive   Robust
      grey  <- c("grey30", "grey70", "grey50")
      #             Boot   Naive   Robust
      lines <- c("dashed","dotted","dotdash") 
      
      
      gCOP_b <- ggplot(dfCO_b, aes(x=Positivity,y=Coverage, col=Type)) + 
        geom_line(aes(linetype=Type),size=1) +
        coord_cartesian(ylim = c(0, 1)) +
        geom_hline(aes(yintercept = 0.95)) +
        xlab("Strenght of practical positivity violation") + 
        scale_x_continuous(breaks=c(0.1, 1, 2),labels=c("Weak", "Moderate", "Strong")) +  
        theme_bw() + theme(legend.position="none", 
                           axis.title = element_text(size = SIZET),
                           axis.text = element_text(size = SIZET),
                           axis.text.x = element_text(hjust = c(0.5,0.5,0.8)) ) +
        ylab("Coverage 95% CI") + 
        scale_color_manual(values=c(grey[1],grey[2],grey[3])) +
        scale_linetype_manual(values = lines) +
        scale_fill_manual(values=c(grey[1],grey[2],grey[3]),
                          name="Method",
                          breaks=unique(dfCO_b$Method),
                          labels=unique(dfCO_b$Method))
      
      gCOP_c <- ggplot(dfCO_c, aes(x=Positivity,y=Coverage, col=Type)) + 
        geom_line(aes(linetype=Type),size=1) +
        coord_cartesian(ylim = c(0, 1)) +
        geom_hline(aes(yintercept = 0.95)) +
        xlab("Strenght of practical positivity violation") + 
        scale_x_continuous(breaks=c(1, 3, 5),labels=c("Weak", "Moderate", "Strong")) + 
        theme_bw() + theme(legend.position="none", 
                           axis.title = element_text(size = SIZET),
                           axis.text = element_text(size = SIZET),
                           axis.text.x = element_text(hjust = c(0.5,0.5,0.8)) ) +
        ylab("Coverage 95% CI") + 
        scale_color_manual(values=c(grey[1],grey[2],grey[3])) +
        scale_linetype_manual(values = lines) +
        scale_fill_manual(values=c(grey[1],grey[2],grey[3]),
                          name="Method",
                          breaks=unique(dfCO_c$Method),
                          labels=unique(dfCO_c$Method))
      
      
      
      
      #########################################################################################################
      #######TIME
      dfTIME_b <- dfTIME[which(dfTIME$Treatment=='binary'),]
      dfTIME_c <- dfTIME[which(dfTIME$Treatment=='continuous'),]
      
      
      #             Boot   Naive   Robust
      grey  <- c("grey30", "grey70", "grey50")
      #             Boot   Naive   Robust
      lines <- c("dashed","dotted","dotdash") 
      
      
      gTIMEP_b <- ggplot(dfTIME_b, aes(x=Positivity,y=Time, col=Type)) + 
        geom_line(aes(linetype=Type),size=1) +
        coord_cartesian(ylim = c(0, 3)) +
        xlab("Strenght of practical positivity violation") + 
        scale_x_continuous(breaks=c(0.1, 1, 2),labels=c("Weak", "Moderate", "Strong")) + 
        theme_bw() + theme(legend.position="none", 
                           axis.title = element_text(size = SIZET),
                           axis.text = element_text(size = SIZET),
                           axis.text.x = element_text(hjust = c(0.5,0.5,0.8)) ) +
        ylab("Time") + 
        scale_color_manual(values=c(grey[1],grey[2],grey[3])) +
        scale_linetype_manual(values = lines) +
        scale_fill_manual(values=c(grey[1],grey[2],grey[3]),
                          name="Method",
                          breaks=unique(dfTIME_b$Method),
                          labels=unique(dfTIME_b$Method))
      
      gTIMEP_c <- ggplot(dfTIME_c, aes(x=Positivity,y=Time, col=Type)) + 
        geom_line(aes(linetype=Type),size=1) +
        coord_cartesian(ylim = c(0, 3)) +
        xlab("Strenght of practical positivity violation") + 
        scale_x_continuous(breaks=c(1, 3, 5),labels=c("Weak", "Moderate", "Strong")) +  
        theme_bw() + theme(legend.position="none", 
                           axis.title = element_text(size = SIZET),
                           axis.text = element_text(size = SIZET),
                           axis.text.x = element_text(hjust = c(0.5,0.5,0.8)) ) +
        ylab("Time") + 
        scale_color_manual(values=c(grey[1],grey[2],grey[3])) +
        scale_linetype_manual(values = lines) +
        scale_fill_manual(values=c(grey[1],grey[2],grey[3]),
                          name="Method",
                          breaks=unique(dfTIME_c$Method),
                          labels=unique(dfTIME_c$Method))
      
      #########################################################################################################
      #########################################################################################################
      # Misspecification
      #########################################################################################################
      #########################################################################################################
      
      
      load(paste("data_results/simu_errors_n",n,"misspecification.Rdata"))
      
      
      #########################################################################################################
      #######SE
      dfSE_b <- dfSE[which(dfSE$Treatment=='binary'),]
      dfSE_c <- dfSE[which(dfSE$Treatment=='continuous'),]
      
      #             Boot  Empirical   Naive   Robust
      grey  <- c("grey30", "grey0", "grey70", "grey50")
      #             Boot  Empirical   Naive   Robust
      lines <- c("dashed","solid","dotted","dotdash") 
      
      dfSE_b$Temp <- dfSE_b$SE[which(dfSE_b$Type=="Empirical")]/dfSE_b$SE
      gSEM_b <- ggplot(dfSE_b, aes(x=Misspecification,y=Temp, col=Type)) + 
        geom_line(aes(linetype=Type),size=1) +
        coord_cartesian(ylim = c(0, 1.2)) +
        xlab("Strenght of misspecification") + 
        scale_x_continuous(breaks=c(0, 0.5, 1),labels=c("Null", "Moderate", "Strong")) + 
        theme_bw() + theme(legend.position="none", 
                           axis.title = element_text(size = SIZET),
                           axis.text = element_text(size = SIZET),
                           axis.text.x = element_text(hjust = c(0.5,0.5,0.8)) ) +
        ylab("SD/SE Ratio") + 
        scale_color_manual(values=c(grey[1],grey[2],grey[3],grey[4])) +
        scale_linetype_manual(values = lines) +
        scale_fill_manual(values=c(grey[1],grey[2],grey[3],grey[4]),
                          name="Method",
                          breaks=unique(dfSE_b$Method),
                          labels=unique(dfSE_b$Method))
      
      dfSE_c$Temp <- dfSE_c$SE[which(dfSE_c$Type=="Empirical")]/dfSE_c$SE
      gSEM_c <- ggplot(dfSE_c, aes(x=Misspecification,y=Temp, col=Type)) + 
        geom_line(aes(linetype=Type),size=1) +
        coord_cartesian(ylim = c(0, 1)) +
        xlab("Strenght of misspecification") + 
        scale_x_continuous(breaks=c(0, 0.5, 1),labels=c("Null", "Moderate", "Strong")) + 
        theme_bw() + theme(legend.position="none", 
                           axis.title = element_text(size = SIZET),
                           axis.text = element_text(size = SIZET),
                           axis.text.x = element_text(hjust = c(0.5,0.5,0.8)) ) +
        ylab("SD/SE Ratio") + 
        scale_color_manual(values=c(grey[1],grey[2],grey[3],grey[4])) +
        scale_linetype_manual(values = lines) +
        scale_fill_manual(values=c(grey[1],grey[2],grey[3],grey[4]),
                          name="Method",
                          breaks=unique(dfSE_c$Method),
                          labels=unique(dfSE_c$Method))
      
      
      #########################################################################################################
      #######COVERAGE
      dfCO_b <- dfCO[which(dfCO$Treatment=='binary'),]
      dfCO_c <- dfCO[which(dfCO$Treatment=='continuous'),]
      
      
      #             Boot   Naive   Robust
      grey  <- c("grey30", "grey70", "grey50")
      #             Boot   Naive   Robust
      lines <- c("dashed","dotted","dotdash") 
      
      
      gCOM_b <- ggplot(dfCO_b, aes(x=Misspecification,y=Coverage, col=Type)) + 
        geom_line(aes(linetype=Type),size=1) +
        coord_cartesian(ylim = c(0, 1)) +
        geom_hline(aes(yintercept = 0.95)) +
        xlab("Strenght of misspecification") + 
        scale_x_continuous(breaks=c(0, 0.5, 1),labels=c("Null", "Moderate", "Strong")) + 
        theme_bw() + theme(legend.position="none", 
                           axis.title = element_text(size = SIZET),
                           axis.text = element_text(size = SIZET),
                           axis.text.x = element_text(hjust = c(0.5,0.5,0.8)) ) +
        ylab("Coverage 95% CI") + 
        scale_color_manual(values=c(grey[1],grey[2],grey[3])) +
        scale_linetype_manual(values = lines) +
        scale_fill_manual(values=c(grey[1],grey[2],grey[3]),
                          name="Method",
                          breaks=unique(dfCO_b$Method),
                          labels=unique(dfCO_b$Method))
      
      gCOM_c <- ggplot(dfCO_c, aes(x=Misspecification,y=Coverage, col=Type)) + 
        geom_line(aes(linetype=Type),size=1) +
        coord_cartesian(ylim = c(0, 1)) +
        geom_hline(aes(yintercept = 0.95)) +
        xlab("Strenght of misspecification") + 
        scale_x_continuous(breaks=c(0, 0.5, 1),labels=c("Null", "Moderate", "Strong")) +
        theme_bw() + theme(legend.position="none", 
                           axis.title = element_text(size = SIZET),
                           axis.text = element_text(size = SIZET),
                           axis.text.x = element_text(hjust = c(0.5,0.5,0.8)) ) +
        ylab("Coverage 95% CI") + 
        scale_color_manual(values=c(grey[1],grey[2],grey[3])) +
        scale_linetype_manual(values = lines) +
        scale_fill_manual(values=c(grey[1],grey[2],grey[3]),
                          name="Method",
                          breaks=unique(dfCO_c$Method),
                          labels=unique(dfCO_c$Method))
      
      
      
      
      #########################################################################################################
      #######TIME
      dfTIME_b <- dfTIME[which(dfTIME$Treatment=='binary'),]
      dfTIME_c <- dfTIME[which(dfTIME$Treatment=='continuous'),]
      
      
      #             Boot   Naive   Robust
      grey  <- c("grey30", "grey70", "grey50")
      #             Boot   Naive   Robust
      lines <- c("dashed","dotted","dotdash") 
      
      
      gTIMEM_b <- ggplot(dfTIME_b, aes(x=Misspecification,y=Time, col=Type)) + 
        geom_line(aes(linetype=Type),size=1) +
        coord_cartesian(ylim = c(0, 3)) +
        xlab("Strenght of misspecification") + 
        scale_x_continuous(breaks=c(0, 0.5, 1),labels=c("Null", "Moderate", "Strong")) +
        theme_bw() + theme(legend.position="none", 
                           axis.title = element_text(size = SIZET),
                           axis.text = element_text(size = SIZET),
                           axis.text.x = element_text(hjust = c(0.5,0.5,0.8)) ) +
        ylab("Time") + 
        scale_color_manual(values=c(grey[1],grey[2],grey[3])) +
        scale_linetype_manual(values = lines) +
        scale_fill_manual(values=c(grey[1],grey[2],grey[3]),
                          name="Method",
                          breaks=unique(dfTIME_b$Method),
                          labels=unique(dfTIME_b$Method))
      
      gTIMEM_c <- ggplot(dfTIME_c, aes(x=Misspecification,y=Time, col=Type)) + 
        geom_line(aes(linetype=Type),size=1) +
        coord_cartesian(ylim = c(0, 3)) +
        xlab("Strenght of misspecification") + 
        scale_x_continuous(breaks=c(0, 0.5, 1),labels=c("Null", "Moderate", "Strong")) +
        theme_bw() + theme(legend.position="none", 
                           axis.title = element_text(size = SIZET),
                           axis.text = element_text(size = SIZET),
                           axis.text.x = element_text(hjust = c(0.5,0.5,0.8)) ) +
        ylab("Time") + 
        scale_color_manual(values=c(grey[1],grey[2],grey[3])) +
        scale_linetype_manual(values = lines) +
        scale_fill_manual(values=c(grey[1],grey[2],grey[3]),
                          name="Method",
                          breaks=unique(dfTIME_c$Method),
                          labels=unique(dfTIME_c$Method))
      
      
      #########################################################################################################
      #########################################################################################################
      # Censoring
      #########################################################################################################
      #########################################################################################################
      
      load(paste("data_results/simu_errors_n",n,"censoring.Rdata"))
      
      
      #########################################################################################################
      #######SE
      dfSE_b <- dfSE[which(dfSE$Treatment=='binary'),]
      dfSE_c <- dfSE[which(dfSE$Treatment=='continuous'),]
      
      #             Boot  Empirical   Naive   Robust
      grey  <- c("grey30", "grey0", "grey70", "grey50")
      #             Boot  Empirical   Naive   Robust
      lines <- c("dashed","solid","dotted","dotdash") 
      
      dfSE_b$Temp <- dfSE_b$SE[which(dfSE_b$Type=="Empirical")]/dfSE_b$SE
      gSEC_b <- ggplot(dfSE_b, aes(x=Censoring,y=Temp, col=Type)) + 
        geom_line(aes(linetype=Type),size=1) +
        coord_cartesian(ylim = c(0, 1)) +
        xlab("Percentage of censoring") + 
        scale_x_continuous(breaks=c(0, 0.4, 0.8),labels=c("0%", "40%", "80%")) +
        theme_bw() + theme(legend.position="none", 
                           axis.title = element_text(size = SIZET),
                           axis.text = element_text(size = SIZET),
                           axis.text.x = element_text(hjust = c(0.5,0.5,0.8)) ) +
        ylab("SD/SE Ratio") + 
        scale_color_manual(values=c(grey[1],grey[2],grey[3],grey[4])) +
        scale_linetype_manual(values = lines) +
        scale_fill_manual(values=c(grey[1],grey[2],grey[3],grey[4]),
                          name="Method",
                          breaks=unique(dfSE_b$Method),
                          labels=unique(dfSE_b$Method))
      
      dfSE_c$Temp <- dfSE_c$SE[which(dfSE_c$Type=="Empirical")]/dfSE_c$SE
      gSEC_c <- ggplot(dfSE_c, aes(x=Censoring,y=Temp, col=Type)) + 
        geom_line(aes(linetype=Type),size=1) +
        coord_cartesian(ylim = c(0, 1)) +
        xlab("Percentage of censoring") + 
        scale_x_continuous(breaks=c(0, 0.4, 0.8),labels=c("0%", "40%", "80%")) +
        theme_bw() + theme(legend.position="none", 
                           axis.title = element_text(size = SIZET),
                           axis.text = element_text(size = SIZET),
                           axis.text.x = element_text(hjust = c(0.5,0.5,0.8)) ) +
        ylab("SD/SE Ratio") + 
        scale_color_manual(values=c(grey[1],grey[2],grey[3],grey[4])) +
        scale_linetype_manual(values = lines) +
        scale_fill_manual(values=c(grey[1],grey[2],grey[3],grey[4]),
                          name="Method",
                          breaks=unique(dfSE_c$Method),
                          labels=unique(dfSE_c$Method))
      
      
      #########################################################################################################
      #######COVERAGE
      dfCO_b <- dfCO[which(dfCO$Treatment=='binary'),]
      dfCO_c <- dfCO[which(dfCO$Treatment=='continuous'),]
      
      
      #             Boot   Naive   Robust
      grey  <- c("grey30", "grey70", "grey50")
      #             Boot   Naive   Robust
      lines <- c("dashed","dotted","dotdash") 
      
      
      gCOC_b <- ggplot(dfCO_b, aes(x=Censoring,y=Coverage, col=Type)) + 
        geom_line(aes(linetype=Type),size=1) +
        coord_cartesian(ylim = c(0, 1)) +
        geom_hline(aes(yintercept = 0.95)) +
        xlab("Percentage of censoring") + 
        scale_x_continuous(breaks=c(0, 0.4, 0.8),labels=c("0%", "40%", "80%")) +
        theme_bw() + theme(legend.position="none", 
                           axis.title = element_text(size = SIZET),
                           axis.text = element_text(size = SIZET),
                           axis.text.x = element_text(hjust = c(0.5,0.5,0.8)) ) +
        ylab("Coverage 95% CI") + 
        scale_color_manual(values=c(grey[1],grey[2],grey[3])) +
        scale_linetype_manual(values = lines) +
        scale_fill_manual(values=c(grey[1],grey[2],grey[3]),
                          name="Method",
                          breaks=unique(dfCO_b$Method),
                          labels=unique(dfCO_b$Method))
      
      gCOC_c <- ggplot(dfCO_c, aes(x=Censoring,y=Coverage, col=Type)) + 
        geom_line(aes(linetype=Type),size=1) +
        coord_cartesian(ylim = c(0, 1)) +
        geom_hline(aes(yintercept = 0.95)) +
        xlab("Censoring") + 
        theme_bw() + theme(legend.position="none", 
                           axis.title = element_text(size = SIZET),
                           axis.text = element_text(size = SIZET),
                           axis.text.x = element_text(hjust = c(0.5,0.5,0.8)) ) +
        ylab("Coverage 95% CI") + 
        scale_color_manual(values=c(grey[1],grey[2],grey[3])) +
        scale_linetype_manual(values = lines) +
        scale_fill_manual(values=c(grey[1],grey[2],grey[3]),
                          name="Method",
                          breaks=unique(dfCO_c$Method),
                          labels=unique(dfCO_c$Method))
      
      
      
      
      #########################################################################################################
      #######TIME
      dfTIME_b <- dfTIME[which(dfTIME$Treatment=='binary'),]
      dfTIME_c <- dfTIME[which(dfTIME$Treatment=='continuous'),]
      
      
      #             Boot   Naive   Robust
      grey  <- c("grey30", "grey70", "grey50")
      #             Boot   Naive   Robust
      lines <- c("dashed","dotted","dotdash") 
      
      
      gTIMEC_b <- ggplot(dfTIME_b, aes(x=Censoring,y=Time, col=Type)) + 
        geom_line(aes(linetype=Type),size=1) +
        coord_cartesian(ylim = c(0, 3)) +
        xlab("Percentage of censoring") + 
        scale_x_continuous(breaks=c(0, 0.4, 0.8),labels=c("0%", "40%", "80%")) +
        theme_bw() + theme(legend.position="none", 
                           axis.title = element_text(size = SIZET),
                           axis.text = element_text(size = SIZET),
                           axis.text.x = element_text(hjust = c(0.5,0.5,0.8)) ) +
        ylab("Time") + 
        scale_color_manual(values=c(grey[1],grey[2],grey[3])) +
        scale_linetype_manual(values = lines) +
        scale_fill_manual(values=c(grey[1],grey[2],grey[3]),
                          name="Method",
                          breaks=unique(dfTIME_b$Method),
                          labels=unique(dfTIME_b$Method))
      
      gTIMEC_c <- ggplot(dfTIME_c, aes(x=Censoring,y=Time, col=Type)) + 
        geom_line(aes(linetype=Type),size=1) +
        coord_cartesian(ylim = c(0, 3)) +
        xlab("Percentage of censoring") + 
        scale_x_continuous(breaks=c(0, 0.4, 0.8),labels=c("0%", "40%", "80%")) + 
        theme_bw() + theme(legend.position="none", 
                           axis.title = element_text(size = SIZET),
                           axis.text = element_text(size = SIZET),
                           axis.text.x = element_text(hjust = c(0.5,0.5,0.8)) ) +
        ylab("Time") + 
        scale_color_manual(values=c(grey[1],grey[2],grey[3])) +
        scale_linetype_manual(values = lines) +
        scale_fill_manual(values=c(grey[1],grey[2],grey[3]),
                          name="Method",
                          breaks=unique(dfTIME_c$Method),
                          labels=unique(dfTIME_c$Method))
      
      
      
      
      
      Noprint <- ggplot(dfTIME_c, aes(x=Censoring,y=Time, col=Type)) + 
        geom_line(aes(linetype=Type),size=1) +
        coord_cartesian(ylim = c(0, 3)) +
        xlab("Percentage of censoring") + 
        scale_x_continuous(breaks=c(0, 0.4, 0.8),labels=c("0%", "40%", "80%")) + 
        theme_bw() + theme(legend.position="bottom", 
                           axis.title = element_text(size = SIZET),
                           axis.text = element_text(size = SIZET),
                           legend.text = element_text(size = SIZELe),
                           legend.title = element_text(size = SIZELe),
                           legend.key.width = unit(LENGTHLe,"cm"),
                           axis.text.x = element_text(hjust = c(0.5,0.5,0.8)) ) +
        ylab("Time") + 
        scale_color_manual(values=c(grey[1],grey[2],grey[3])) +
        scale_linetype_manual(values = lines) +
        scale_fill_manual(values=c(grey[1],grey[2],grey[3]),
                          name="Method",
                          breaks=unique(dfTIME_c$Method),
                          labels=unique(dfTIME_c$Method))
      
      
      
      
      #########################################################################################################
      #######BINARY TRT - Bias, MSE across positivity, misspecification and censoring
      # Extract the legend. Returns a gtable
      leg <- get_legend(Noprint)
      
      
      grid.newpage()
      pushViewport(viewport(layout = grid.layout(10, 2, heights = unit(c(0.5, 0.5, 5, 
                                                                         0.5, 0.5, 5, 
                                                                         0.5, 0.5, 5, 
                                                                         1), "null"))))  
      
      grid.text("Practical positivity violation", 
                vp = viewport(layout.pos.row = 1, layout.pos.col = 1:2),gp=gpar(fontsize=20))
      
      
      print(gSEP_b, vp = viewport(layout.pos.row = 3, layout.pos.col = 1))
      print(gCOP_b, vp = viewport(layout.pos.row = 3, layout.pos.col = 2))
      
      grid.text("Misspecification", vp = viewport(layout.pos.row = 4, layout.pos.col = 1:2),gp=gpar(fontsize=20))
      
      print(gSEM_b, vp = viewport(layout.pos.row = 6, layout.pos.col = 1))
      print(gCOM_b, vp = viewport(layout.pos.row = 6, layout.pos.col = 2))
      
      grid.text("Censoring", vp = viewport(layout.pos.row = 7, layout.pos.col = 1:2),gp=gpar(fontsize=20))
      
      print(gSEC_b, vp = viewport(layout.pos.row = 9, layout.pos.col = 1))
      print(gCOC_b, vp = viewport(layout.pos.row = 9, layout.pos.col = 2))
      
      print(as_ggplot(leg), vp = viewport(layout.pos.row = 10, layout.pos.col = 1:2),gp=gpar(fontsize=20))

      
      
      
      
      #########################################################################################################
      #######CONTINUOUS TRT - Bias, MSE across positivity, misspecification and censoring
      # Extract the legend. Returns a gtable
      leg <- get_legend(Noprint_c)
      
      
      grid.newpage()
      pushViewport(viewport(layout = grid.layout(10, 2, heights = unit(c(0.5, 0.5, 5, 
                                                                         0.5, 0.5, 5, 
                                                                         0.5, 0.5, 5, 
                                                                         1), "null"))))  
      
      grid.text("Practical positivity violation", 
                vp = viewport(layout.pos.row = 1, layout.pos.col = 1:2),gp=gpar(fontsize=20))
      
      
      print(gSEP_c, vp = viewport(layout.pos.row = 3, layout.pos.col = 1))
      print(gCOP_c, vp = viewport(layout.pos.row = 3, layout.pos.col = 2))
      
      grid.text("Misspecification", vp = viewport(layout.pos.row = 4, layout.pos.col = 1:2),gp=gpar(fontsize=20))
      
      print(gSEM_c, vp = viewport(layout.pos.row = 6, layout.pos.col = 1))
      print(gCOM_c, vp = viewport(layout.pos.row = 6, layout.pos.col = 2))
      
      grid.text("Censoring", vp = viewport(layout.pos.row = 7, layout.pos.col = 1:2),gp=gpar(fontsize=20))
      
      print(gSEC_c, vp = viewport(layout.pos.row = 9, layout.pos.col = 1))
      print(gCOC_c, vp = viewport(layout.pos.row = 9, layout.pos.col = 2))
      
      print(as_ggplot(leg), vp = viewport(layout.pos.row = 10, layout.pos.col = 1:2))





###################################################################################################################################################
###################################################################################################################################################
#
#
# Simulation OVER N and OVER COVARIATES
#
#
###################################################################################################################################################
###################################################################################################################################################

      SIZEL <- 2
      SIZET <- 14
      SIZEP <- 4

      #########################################################################################################
      #########################################################################################################
      # Sample Size - N
      #########################################################################################################
      #########################################################################################################
      
      load(paste("data_results/simu_n_sample_size.Rdata"))
      
      #########################################################################################################
      #######BIAS
      dfBIAS_b <- dfBIAS[which(dfBIAS$Treatment=='binary'),]
      dfBIAS_c <- dfBIAS[which(dfBIAS$Treatment=='continuous'),]
      
      gBiasN_b <- ggplot(dfBIAS_b, aes(x=`Sample size`,y=abs(Bias), col=Method)) + 
        geom_line(aes(linetype=Method),size=SIZEL) +
        geom_point(aes(linetype=Method),size=SIZEP) +
        xlab("Sample Size") + 
        coord_cartesian(ylim = c(0, 2)) +
        scale_color_manual(values=c("black")) +
        theme_bw() + theme(legend.position="none", 
                           axis.title = element_text(size = SIZET),
                           axis.text = element_text(size = SIZET),
                           axis.text.x = element_text(hjust = c(0.5,0.5,0.8)) ) +
        ylab("Abs(Bias)")  
      
      
      gBiasN_c <- ggplot(dfBIAS_c, aes(x=`Sample size`,y=abs(Bias), col=Method)) + 
        geom_line(aes(linetype=Method),size=SIZEL) +
        geom_point(aes(linetype=Method),size=SIZEP) +
        xlab("Sample Size") + 
        coord_cartesian(ylim = c(0, 1)) +
        scale_color_manual(values=c("black")) +
        theme_bw() + theme(legend.position="none", 
                           axis.title = element_text(size = SIZET),
                           axis.text = element_text(size = SIZET),
                           axis.text.x = element_text(hjust = c(0.5,0.5,0.8)) ) +
        ylab("Abs(Bias)")  
      
      
      #########################################################################################################
      #######MSE
      dfMSE_b <- dfMSE[which(dfMSE$Treatment=='binary'),]
      dfMSE_c <- dfMSE[which(dfMSE$Treatment=='continuous'),]
      
      
      gMseN_b <- ggplot(dfMSE_b, aes(x=`Sample size`,y=sqrt(MSE), col=Method)) + 
        geom_line(aes(linetype=Method),size=SIZEL) +
        geom_point(aes(linetype=Method),size=SIZEP) +
        xlab("Sample Size") + 
        coord_cartesian(ylim = c(0, 2)) +
        scale_color_manual(values=c("black")) +
        theme_bw() + theme(legend.position="none", 
                           axis.title = element_text(size = SIZET),
                           axis.text = element_text(size = SIZET),
                           axis.text.x = element_text(hjust = c(0.5,0.5,0.8)) ) +
        ylab("RMSE")  
      
      
      gMseN_c <- ggplot(dfMSE_c, aes(x=`Sample size`,y=sqrt(MSE), col=Method)) + 
        geom_line(aes(linetype=Method),size=SIZEL) +
        geom_point(aes(linetype=Method),size=SIZEP) +
        xlab("Sample Size") + 
        coord_cartesian(ylim = c(0, 2)) +
        scale_color_manual(values=c("black")) +
        theme_bw() + theme(legend.position="none", 
                           axis.title = element_text(size = SIZET),
                           axis.text = element_text(size = SIZET),
                           axis.text.x = element_text(hjust = c(0.5,0.5,0.8)) ) +
        ylab("RMSE")  
      
      
      
      #########################################################################################################
      #######BAL
      dfBAL_b <- dfBAL[which(dfBAL$Treatment=='binary'),]
      dfBAL_c <- dfBAL[which(dfBAL$Treatment=='continuous'),]
      
      
      gBalN_b <- ggplot(dfBAL_b, aes(x=`Sample size`,y=Balance, col=Method)) + 
        geom_line(aes(linetype=Method),size=SIZEL) +
        geom_point(aes(linetype=Method),size=SIZEP) +
        xlab("Sample Size") + 
        coord_cartesian(ylim = c(0, 0.1)) +
        scale_color_manual(values=c("black")) +
        theme_bw() + theme(legend.position="none", 
                           axis.title = element_text(size = SIZET),
                           axis.text = element_text(size = SIZET),
                           axis.text.x = element_text(hjust = c(0.5,0.5,0.8)) ) +
        ylab("Balance")  
      
      
      gBalN_c <- ggplot(dfBAL_c, aes(x=`Sample size`,y=Balance, col=Method)) + 
        geom_line(aes(linetype=Method),size=SIZEL) +
        geom_point(aes(linetype=Method),size=SIZEP) +
        xlab("Sample Size") + 
        coord_cartesian(ylim = c(0, 0.1)) +
        scale_color_manual(values=c("black")) +
        theme_bw() + theme(legend.position="none", 
                           axis.title = element_text(size = SIZET),
                           axis.text = element_text(size = SIZET),
                           axis.text.x = element_text(hjust = c(0.5,0.5,0.8)) ) +
        ylab("Balance")  
      
      
      
      #########################################################################################################
      #######TIME
      dfTIME_b <- dfTIME[which(dfTIME$Treatment=='binary'),]
      dfTIME_c <- dfTIME[which(dfTIME$Treatment=='continuous'),]
      
      
      gTimeN_b <- ggplot(dfTIME_b, aes(x=`Sample size`,y=Time, col=Method)) + 
        geom_line(aes(linetype=Method),size=SIZEL) +
        geom_point(aes(linetype=Method),size=SIZEP) +
        xlab("Sample Size") + 
        coord_cartesian(ylim = c(0, 1.5)) +
        scale_color_manual(values=c("black")) +
        theme_bw() + theme(legend.position="none", 
                           axis.title = element_text(size = SIZET),
                           axis.text = element_text(size = SIZET),
                           axis.text.x = element_text(hjust = c(0.5,0.5,0.8)) ) +
        ylab("Average time in secs")  
      
      
      gTimeN_c <- ggplot(dfTIME_c, aes(x=`Sample size`,y=Time, col=Method)) + 
        geom_line(aes(linetype=Method),size=SIZEL) +
        geom_point(aes(linetype=Method),size=SIZEP) +
        xlab("Sample Size") + 
        coord_cartesian(ylim = c(0, 1.5)) +
        scale_color_manual(values=c("black")) +
        theme_bw() + theme(legend.position="none", 
                           axis.title = element_text(size = SIZET),
                           axis.text = element_text(size = SIZET),
                           axis.text.x = element_text(hjust = c(0.5,0.5,0.8)) ) +
        ylab("Average time in secs")  
      
      
      
      
      #########################################################################################################
      #######BINARY TRT - Bias, MSE across positivity, misspecification and censoring
      # Extract the legend. Returns a gtable

      
      grid.newpage()
      pushViewport(viewport(layout = grid.layout(7, 2, heights = unit(c(0.5, 0.5, 5, 
                                                                         0.5, 0.5, 5, 
                                                                         1), "null"))))  
      
      grid.text("Binary treatment", 
                vp = viewport(layout.pos.row = 1, layout.pos.col = 1:2),gp=gpar(fontsize=20))
      
      
      print(gBiasN_b, vp = viewport(layout.pos.row = 3, layout.pos.col = 1))
      print(gMseN_b, vp = viewport(layout.pos.row = 3, layout.pos.col = 2))
      
      print(gBalN_b, vp = viewport(layout.pos.row = 6, layout.pos.col = 1))
      print(gTimeN_b, vp = viewport(layout.pos.row = 6, layout.pos.col = 2))
      
      
      
      #########################################################################################################
      #######CONTINUOUS TRT - Bias, MSE across positivity, misspecification and censoring
      # Extract the legend. Returns a gtable
      
      
      grid.newpage()
      pushViewport(viewport(layout = grid.layout(7, 2, heights = unit(c(0.5, 0.5, 5, 
                                                                        0.5, 0.5, 5, 
                                                                        1), "null"))))  
      
      grid.text("Continuous treatment", 
                vp = viewport(layout.pos.row = 1, layout.pos.col = 1:2),gp=gpar(fontsize=20))
      
      
      print(gBiasN_c, vp = viewport(layout.pos.row = 3, layout.pos.col = 1))
      print(gMseN_c, vp = viewport(layout.pos.row = 3, layout.pos.col = 2))
      
      print(gBalN_c, vp = viewport(layout.pos.row = 6, layout.pos.col = 1))
      print(gTimeN_c, vp = viewport(layout.pos.row = 6, layout.pos.col = 2))
      
      
      
      
      
      
      
      
      #########################################################################################################
      #########################################################################################################
      # Number Covariates
      #########################################################################################################
      #########################################################################################################
      
      SIZEL <- 2
      SIZET <- 14
      SIZEP <- 4
      
      n <- 1000
      
      load(paste("data_results/simu_n",n,"num_cova.Rdata"))
      
      #########################################################################################################
      #######BIAS
      dfBIAS_b <- dfBIAS[which(dfBIAS$Treatment=='binary'),]
      dfBIAS_c <- dfBIAS[which(dfBIAS$Treatment=='continuous'),]
      
      gBiasC_b <- ggplot(dfBIAS_b, aes(x=`Number Covariates`,y=abs(Bias), col=Method)) + 
        geom_line(aes(linetype=Method),size=SIZEL) +
        geom_point(aes(linetype=Method),size=SIZEP) +
        xlab("Number covariates") + 
        coord_cartesian(ylim = c(0, 0.2)) +
        scale_color_manual(values=c("black")) +
        theme_bw() + theme(legend.position="none", 
                           axis.title = element_text(size = SIZET),
                           axis.text = element_text(size = SIZET),
                           axis.text.x = element_text(hjust = c(0.5,0.5,0.8)) ) +
        ylab("Abs(Bias)")  
      
      
      gBiasC_c <- ggplot(dfBIAS_c, aes(x=`Number Covariates`,y=abs(Bias), col=Method)) + 
        geom_line(aes(linetype=Method),size=SIZEL) +
        geom_point(aes(linetype=Method),size=SIZEP) +
        xlab("Number covariates") + 
        coord_cartesian(ylim = c(0, 0.1)) +
        scale_color_manual(values=c("black")) +
        theme_bw() + theme(legend.position="none", 
                           axis.title = element_text(size = SIZET),
                           axis.text = element_text(size = SIZET),
                           axis.text.x = element_text(hjust = c(0.5,0.5,0.8)) ) +
        ylab("Abs(Bias)")  
      
      
      
      #########################################################################################################
      #######MSE
      dfMSE_b <- dfMSE[which(dfMSE$Treatment=='binary'),]
      dfMSE_c <- dfMSE[which(dfMSE$Treatment=='continuous'),]
      
      
      gMseC_b <- ggplot(dfMSE_b, aes(x=`Number Covariates`,y=sqrt(MSE), col=Method)) + 
        geom_line(aes(linetype=Method),size=SIZEL) +
        geom_point(aes(linetype=Method),size=SIZEP) +
        xlab("Number covariates") + 
        coord_cartesian(ylim = c(0, 0.25)) +
        scale_color_manual(values=c("black")) +
        theme_bw() + theme(legend.position="none", 
                           axis.title = element_text(size = SIZET),
                           axis.text = element_text(size = SIZET),
                           axis.text.x = element_text(hjust = c(0.5,0.5,0.8)) ) +
        ylab("RMSE")  
      
      
      gMseC_c <- ggplot(dfMSE_c, aes(x=`Number Covariates`,y=sqrt(MSE), col=Method)) + 
        geom_line(aes(linetype=Method),size=SIZEL) +
        geom_point(aes(linetype=Method),size=SIZEP) +
        xlab("Number covariates") + 
        coord_cartesian(ylim = c(0, 0.1)) +
        scale_color_manual(values=c("black")) +
        theme_bw() + theme(legend.position="none", 
                           axis.title = element_text(size = SIZET),
                           axis.text = element_text(size = SIZET),
                           axis.text.x = element_text(hjust = c(0.5,0.5,0.8)) ) +
        ylab("RMSE")  
      
      
      
      #########################################################################################################
      #######BAL
      dfBAL_b <- dfBAL[which(dfBAL$Treatment=='binary'),]
      dfBAL_c <- dfBAL[which(dfBAL$Treatment=='continuous'),]
      
      
      gBalC_b <- ggplot(dfBAL_b, aes(x=`Number Covariates`,y=Balance, col=Method)) + 
        geom_line(aes(linetype=Method),size=SIZEL) +
        geom_point(aes(linetype=Method),size=SIZEP) +
        xlab("Number covariates") + 
        coord_cartesian(ylim = c(0, 1)) +
        scale_color_manual(values=c("black")) +
        theme_bw() + theme(legend.position="none", 
                           axis.title = element_text(size = SIZET),
                           axis.text = element_text(size = SIZET),
                           axis.text.x = element_text(hjust = c(0.5,0.5,0.8)) ) +
        ylab("Balance")  
      
      
      gBalC_c <- ggplot(dfBAL_c, aes(x=`Number Covariates`,y=Balance, col=Method)) + 
        geom_line(aes(linetype=Method),size=SIZEL) +
        geom_point(aes(linetype=Method),size=SIZEP) +
        xlab("Number covariates") + 
        coord_cartesian(ylim = c(0, 1)) +
        scale_color_manual(values=c("black")) +
        theme_bw() + theme(legend.position="none", 
                           axis.title = element_text(size = SIZET),
                           axis.text = element_text(size = SIZET),
                           axis.text.x = element_text(hjust = c(0.5,0.5,0.8)) ) +
        ylab("Balance")  
      
      
      
      #########################################################################################################
      #######TIME
      dfTIME_b <- dfTIME[which(dfTIME$Treatment=='binary'),]
      dfTIME_c <- dfTIME[which(dfTIME$Treatment=='continuous'),]
      
      
      gTimeC_b <- ggplot(dfTIME_b, aes(x=`Number Covariates`,y=Time, col=Method)) + 
        geom_line(aes(linetype=Method),size=SIZEL) +
        geom_point(aes(linetype=Method),size=SIZEP) +
        xlab("Number covariates") + 
        coord_cartesian(ylim = c(0, 0.25)) +
        scale_color_manual(values=c("black")) +
        theme_bw() + theme(legend.position="none", 
                           axis.title = element_text(size = SIZET),
                           axis.text = element_text(size = SIZET),
                           axis.text.x = element_text(hjust = c(0.5,0.5,0.8)) ) +
        ylab("Average time in secs")  
      
      
      gTimeC_c <- ggplot(dfTIME_c, aes(x=`Number Covariates`,y=Time, col=Method)) + 
        geom_line(aes(linetype=Method),size=SIZEL) +
        geom_point(aes(linetype=Method),size=SIZEP) +
        xlab("Number covariates") + 
        coord_cartesian(ylim = c(0, 0.25)) +
        scale_color_manual(values=c("black")) +
        theme_bw() + theme(legend.position="none", 
                           axis.title = element_text(size = SIZET),
                           axis.text = element_text(size = SIZET),
                           axis.text.x = element_text(hjust = c(0.5,0.5,0.8)) ) +
        ylab("Average time in secs")  
      
      
      
      
      
      
      #########################################################################################################
      #######BINARY TRT - Bias, MSE across positivity, misspecification and censoring
      # Extract the legend. Returns a gtable
      
      
      grid.newpage()
      pushViewport(viewport(layout = grid.layout(7, 2, heights = unit(c(0.5, 0.5, 5, 
                                                                        0.5, 0.5, 5, 
                                                                        1), "null"))))  
      
      grid.text("Binary treatment", 
                vp = viewport(layout.pos.row = 1, layout.pos.col = 1:2),gp=gpar(fontsize=20))
      
      
      print(gBiasC_b, vp = viewport(layout.pos.row = 3, layout.pos.col = 1))
      print(gMseC_b, vp = viewport(layout.pos.row = 3, layout.pos.col = 2))
      
      print(gBalC_b, vp = viewport(layout.pos.row = 6, layout.pos.col = 1))
      print(gTimeC_b, vp = viewport(layout.pos.row = 6, layout.pos.col = 2))
      
      
      
      #########################################################################################################
      #######CONTINUOUS TRT - Bias, MSE across positivity, misspecification and censoring
      # Extract the legend. Returns a gtable
      
      
      grid.newpage()
      pushViewport(viewport(layout = grid.layout(7, 2, heights = unit(c(0.5, 0.5, 5, 
                                                                        0.5, 0.5, 5, 
                                                                        1), "null"))))  
      
      grid.text("Continuous treatment", 
                vp = viewport(layout.pos.row = 1, layout.pos.col = 1:2),gp=gpar(fontsize=20))
      
      
      print(gBiasC_c, vp = viewport(layout.pos.row = 3, layout.pos.col = 1))
      print(gMseC_c, vp = viewport(layout.pos.row = 3, layout.pos.col = 2))
      
      print(gBalC_c, vp = viewport(layout.pos.row = 6, layout.pos.col = 1))
      print(gTimeC_c, vp = viewport(layout.pos.row = 6, layout.pos.col = 2))
      
      
      
      

###################################################################################################################################################
###################################################################################################################################################
#
#
# Simulation OVER DELTA
#
#
###################################################################################################################################################
###################################################################################################################################################

      
      
      #########################################################################################################
      #########################################################################################################
      # Delta
      #########################################################################################################
      #########################################################################################################
      n <- 1000
      
      load(paste("data_results/simu_n",n,"delta.Rdata"))
      
      
      #########################################################################################################
      #######BIAS
      dfBIAS_b <- dfBIAS[which(dfBIAS$Treatment=='binary' & dfBIAS$Delta<=0.1) ,]
      dfBIAS_c <- dfBIAS[which(dfBIAS$Treatment=='continuous' & dfBIAS$Delta<=0.1),]
      
      gBiasD_b <- ggplot(dfBIAS_b, aes(x=(Delta),y=(abs(Bias)), col=Method)) + 
        geom_line(aes(linetype=Method),size=1) +
        geom_point(aes(linetype=Method),size=3) +
        xlab(expression(paste(delta))) + 
        scale_x_continuous(breaks=c(0, 0.05, 0.1),labels=c("0.0001", "0.05", "0.1")) +
        coord_cartesian(ylim = c(0, 0.8)) +
        scale_color_manual(values=c("black")) +
        theme_bw() + theme(legend.position="none") +
        ylab("Abs(Bias)")  
      
      
      gBiasD_c <- ggplot(dfBIAS_c, aes(x=Delta,y=abs(Bias), col=Method)) + 
        geom_line(aes(linetype=Method),size=1) +
        geom_point(aes(linetype=Method),size=3) +
        xlab(expression(paste(delta))) + 
        scale_x_continuous(breaks=c(0, 0.05, 0.1),labels=c("0.0001", "0.05", "0.1")) +
        coord_cartesian(ylim = c(0, 0.6)) +
        scale_color_manual(values=c("black")) +
        theme_bw() + theme(legend.position="none") +
        ylab("Abs(Bias)") 
      
      
      
      #########################################################################################################
      #######MSE
      dfMSE_b <- dfMSE[which(dfMSE$Treatment=='binary' & dfBIAS$Delta<=0.1),]
      dfMSE_c <- dfMSE[which(dfMSE$Treatment=='continuous' & dfBIAS$Delta<=0.1),]
      
      
      gMseD_b <- ggplot(dfMSE_b, aes(x=Delta,y=sqrt(MSE), col=Method)) + 
        geom_line(aes(linetype=Method),size=1) +
        geom_point(aes(linetype=Method),size=3) +
        xlab(expression(paste(delta))) + 
        scale_x_continuous(breaks=c(0, 0.05, 0.1),labels=c("0.0001", "0.05", "0.1")) +
        coord_cartesian(ylim = c(0, 0.8)) +
        scale_color_manual(values=c("black")) +
        theme_bw() + theme(legend.position="none") +
        ylab("RMSE")  
      
      
      gMseD_c <- ggplot(dfMSE_c, aes(x=Delta,y=sqrt(MSE), col=Method)) + 
        geom_line(aes(linetype=Method),size=1) +
        geom_point(aes(linetype=Method),size=3) +
        xlab(expression(paste(delta))) + 
        scale_x_continuous(breaks=c(0, 0.05, 0.1),labels=c("0.0001", "0.05", "0.1")) +
        coord_cartesian(ylim = c(0, 0.6)) +
        scale_color_manual(values=c("black")) +
        theme_bw() + theme(legend.position="none") +
        ylab("RMSE")  
      
      
      
      #########################################################################################################
      #######BAL
      dfBAL_b <- dfBAL[which(dfBAL$Treatment=='binary' & dfBIAS$Delta<=0.1),]
      dfBAL_c <- dfBAL[which(dfBAL$Treatment=='continuous' & dfBIAS$Delta<=0.1),]
      
      
      gBalD_b <- ggplot(dfBAL_b, aes(x=Delta,y=Balance, col=Method)) + 
        geom_line(aes(linetype=Method),size=1) +
        geom_point(aes(linetype=Method),size=3) +
        xlab(expression(paste(delta))) + 
        scale_x_continuous(breaks=c(0, 0.05, 0.1),labels=c("0.0001", "0.05", "0.1")) +
        coord_cartesian(ylim = c(0, 1)) +
        scale_color_manual(values=c("black")) +
        theme_bw() + theme(legend.position="none") +
        ylab("Balance")  
      
      
      gBalD_c <- ggplot(dfBAL_c, aes(x=Delta,y=Balance, col=Method)) + 
        geom_line(aes(linetype=Method),size=1) +
        geom_point(aes(linetype=Method),size=3) +
        xlab(expression(paste(delta))) +         
        scale_x_continuous(breaks=c(0, 0.05, 0.1),labels=c("0.0001", "0.05", "0.1")) +
        coord_cartesian(ylim = c(0, 1)) +
        scale_color_manual(values=c("black")) +
        theme_bw() + theme(legend.position="none") +
        ylab("Balance")  
      
      
      
      #########################################################################################################
      #######TIME
      dfTIME_b <- dfTIME[which(dfTIME$Treatment=='binary'& dfBIAS$Delta<=0.1),]
      dfTIME_c <- dfTIME[which(dfTIME$Treatment=='continuous'& dfBIAS$Delta<=0.1),]
      
      
      gTimeD_b <- ggplot(dfTIME_b, aes(x=Delta,y=Time, col=Method)) + 
        geom_line(aes(linetype=Method),size=1) +
        geom_point(aes(linetype=Method),size=3) +
        xlab(expression(paste(delta))) + 
        scale_x_continuous(breaks=c(0, 0.05, 0.1),labels=c("0.0001", "0.05", "0.1")) +
        coord_cartesian(ylim = c(0, 0.1)) +
        scale_color_manual(values=c("black")) +
        theme_bw() + theme(legend.position="none") +
        ylab("Time")  
      
      
      gTimeD_c <- ggplot(dfTIME_c, aes(x=Delta,y=Time, col=Method)) + 
        geom_line(aes(linetype=Method),size=1) +
        geom_point(aes(linetype=Method),size=3) +
        xlab(expression(paste(delta))) + 
        scale_x_continuous(breaks=c(0, 0.05, 0.1),labels=c("0.0001", "0.05", "0.1")) +
        coord_cartesian(ylim = c(0, 0.1)) +
        scale_color_manual(values=c("black")) +
        theme_bw() + theme(legend.position="none") +
        ylab("Time") 
      
      
      
      
      #########################################################################################################
      #######BINARY TRT - Bias, MSE across positivity, misspecification and censoring
      # Extract the legend. Returns a gtable
      
      
      grid.newpage()
      pushViewport(viewport(layout = grid.layout(9, 3, heights = unit(c(0.5, 0.5, 0.5, 5, 
                                                                        0.5, 0.5, 0.5, 5, 
                                                                        1), "null"))))  
      grid.text("Binary treatment", vp = viewport(layout.pos.row = 2, layout.pos.col = 1:3))
      
      
      print(gBiasD_b, vp = viewport(layout.pos.row = 4, layout.pos.col = 1))
      print(gMseD_b, vp = viewport(layout.pos.row = 4, layout.pos.col = 2))
      print(gBalD_b, vp = viewport(layout.pos.row = 4, layout.pos.col = 3))
      
      grid.text("Continuous treatment", vp = viewport(layout.pos.row = 6, layout.pos.col = 1:3))
      
      print(gBiasD_c, vp = viewport(layout.pos.row = 8, layout.pos.col = 1))
      print(gMseD_c, vp = viewport(layout.pos.row = 8, layout.pos.col = 2))
      print(gBalD_c, vp = viewport(layout.pos.row = 8, layout.pos.col = 3))
      

    
      
      
      
      #########################################################################################################
      #######CONTINUOUS TRT - Bias, MSE across positivity, misspecification and censoring
      # Extract the legend. Returns a gtable
      
      
      grid.newpage()
      pushViewport(viewport(layout = grid.layout(7, 2, heights = unit(c(0.5, 0.5, 5, 
                                                                        0.5, 0.5, 5, 
                                                                        1), "null"))))  
      
      grid.text("Sample Size", 
                vp = viewport(layout.pos.row = 1, layout.pos.col = 1:2))
      
      
      print(gBiasD_c, vp = viewport(layout.pos.row = 3, layout.pos.col = 1))
      print(gMseD_c, vp = viewport(layout.pos.row = 3, layout.pos.col = 2))
      
      print(gBalD_c, vp = viewport(layout.pos.row = 6, layout.pos.col = 1))
      print(gTimeD_c, vp = viewport(layout.pos.row = 6, layout.pos.col = 2))

