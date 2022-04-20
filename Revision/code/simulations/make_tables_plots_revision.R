########################################################################
########################################################################
#***********************************************************************
#
#   In this file, we provide code to make tables in the revised
#   manuscript and in the supplementary material.
#
#***********************************************************************
########################################################################
########################################################################


library(ggplot2)
library(ggpubr)
library(dplyr)
library(grid)

rm(list = ls())


n <- 1000
D <- 3


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


setwd("/Users/santam13/Documents/NYU/ROW/ROW-time-to-event/code/simulations/")


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

  load(paste("data_results/simu_n",n,"positivity.Rdata"))

  #########################################################################################################
  #######BIAS
  dfBIAS_b <- dfBIAS[which(dfBIAS$Treatment=='binary'),]
  dfBIAS_b$id <- sort(rep(1:(nrow(dfBIAS_b)/2),2))
  dfBIAS_b <- reshape(dfBIAS_b,timevar ="Positivity",direction="wide",idvar="id")
  
  dfBIAS_c <- dfBIAS[which(dfBIAS$Treatment=='continuous'),]
  dfBIAS_c$id <- sort(rep(1:(nrow(dfBIAS_c)/2),2))
  dfBIAS_c <- reshape(dfBIAS_c,timevar ="Positivity",direction="wide",idvar="id")
  
  #########################################################################################################
  #######MSE
  dfMSE_b <- dfMSE[which(dfMSE$Treatment=='binary'),]
  dfMSE_b$id <- sort(rep(1:(nrow(dfMSE_b)/2),2))
  dfMSE_b <- reshape(dfMSE_b,timevar ="Positivity",direction="wide",idvar="id")
  
  dfMSE_c <- dfMSE[which(dfMSE$Treatment=='continuous'),]
  dfMSE_c$id <- sort(rep(1:(nrow(dfMSE_c)/2),2))
  dfMSE_c <- reshape(dfMSE_c,timevar ="Positivity",direction="wide",idvar="id")
  
  #########################################################################################################
  #######TIME
  dfTIME_b <- dfTIME[which(dfTIME$Treatment=='binary'),]
  dfTIME_b$id <- sort(rep(1:(nrow(dfTIME_b)/2),2))
  dfTIME_b <- reshape(dfTIME_b,timevar ="Positivity",direction="wide",idvar="id")
  
  dfTIME_c <- dfTIME[which(dfTIME$Treatment=='continuous'),]
  dfTIME_c$id <- sort(rep(1:(nrow(dfTIME_c)/2),2))
  dfTIME_c <- reshape(dfTIME_c,timevar ="Positivity",direction="wide",idvar="id")
  
  #########################################################################################################
  #######BALANCE
  dfBAL_b <- dfBAL[which(dfBAL$Treatment=='binary'),]
  dfBAL_b$id <- sort(rep(1:(nrow(dfBAL_b)/2),2))
  dfBAL_b <- reshape(dfBAL_b,timevar ="Positivity",direction="wide",idvar="id")
  
  dfBAL_c <- dfBAL[which(dfBAL$Treatment=='continuous'),]
  dfBAL_c$id <- sort(rep(1:(nrow(dfBAL_c)/2),2))
  dfBAL_c <- reshape(dfBAL_c,timevar ="Positivity",direction="wide",idvar="id")
  
  
  #########################################################################################################
  #######Binary treatment - Table
  temp_table_b <- data.frame(dfBIAS_b$Method.1,
                           round(dfBIAS_b$`Abs-Bias.1`,D),round(dfMSE_b$RMSE.1,D), #Correct - no PPV
                           round(dfBIAS_b$`Abs-Bias.2`,D),round(dfMSE_b$RMSE.2,D)) #Correct - PPV 
  
  
  temp2_table_b <- data.frame(dfBIAS_b$Method.1,
                             round(dfBAL_b$Balance.1,D),round(dfTIME_b$Time.1,D), #Correct - no PPV
                             round(dfBAL_b$Balance.2,D),round(dfTIME_b$Time.2,D)) #Correct - PPV 
  
  
  #########################################################################################################
  #######Continuous treatment - Table
  temp_table_c <- data.frame(dfBIAS_c$Method.1,
                             round(dfBIAS_c$`Abs-Bias.1`,D),round(dfMSE_c$RMSE.1,D), #Correct - no PPV
                             round(dfBIAS_c$`Abs-Bias.2`,D),round(dfMSE_c$RMSE.2,D)) #Correct - PPV 
  
  temp2_table_c <- data.frame(dfBIAS_c$Method.1,
                              round(dfBAL_c$Balance.1,D),round(dfTIME_c$Time.1,D), #Correct - no PPV
                              round(dfBAL_c$Balance.2,D),round(dfTIME_c$Time.2,D)) #Correct - PPV 
  
  
  #########################################################################################################
  #########################################################################################################
  # Misspecification
  #########################################################################################################
  #########################################################################################################
  
  load(paste("data_results/simu_n",n,"misspecified.Rdata"))
  
  #########################################################################################################
  #######BIAS
  dfBIAS_b <- dfBIAS[which(dfBIAS$Treatment=='binary'),]
  dfBIAS_b$id <- rep(1:2,(nrow(dfBIAS_b)/2))
  dfBIAS_b <- dfBIAS_b %>% filter(id==1)
  
  dfBIAS_c <- dfBIAS[which(dfBIAS$Treatment=='continuous'),]
  dfBIAS_c$id <- rep(1:2,(nrow(dfBIAS_c)/2))
  dfBIAS_c <- dfBIAS_c %>% filter(id==1)
  
  #########################################################################################################
  #######MSE
  dfMSE_b <- dfMSE[which(dfMSE$Treatment=='binary'),]
  dfMSE_b$id <- rep(1:2,(nrow(dfMSE_b)/2))
  dfMSE_b <- dfMSE_b %>% filter(id==1)
  
  dfMSE_c <- dfMSE[which(dfMSE$Treatment=='continuous'),]
  dfMSE_c$id <- rep(1:2,(nrow(dfMSE_c)/2))
  dfMSE_c <- dfMSE_c %>% filter(id==1)
  
  #########################################################################################################
  #######TIME
  dfTIME_b <- dfTIME[which(dfTIME$Treatment=='binary'),]
  dfTIME_b$id <- rep(1:2,(nrow(dfTIME_b)/2))
  dfTIME_b <- dfTIME_b %>% filter(id==1)
  
  dfTIME_c <- dfTIME[which(dfTIME$Treatment=='continuous'),]
  dfTIME_c$id <- rep(1:2,(nrow(dfTIME_c)/2))
  dfTIME_c <- dfTIME_c %>% filter(id==1)
  
  #########################################################################################################
  #######BALANCE
  dfBAL_b <- dfBAL[which(dfBAL$Treatment=='binary'),]
  dfBAL_b$id <- rep(1:2,(nrow(dfBAL_b)/2))
  dfBAL_b <- dfBAL_b %>% filter(id==1)
  
  dfBAL_c <- dfBAL[which(dfBAL$Treatment=='continuous'),]
  dfBAL_c$id <- rep(1:2,(nrow(dfBAL_c)/2))
  dfBAL_c <- dfBAL_c %>% filter(id==1)
  
  

  #########################################################################################################
  #######Binary treatment - Table
  binary_t_table <- data.frame(temp_table_b,
                               round(dfBIAS_b$`Abs-Bias`,D),round(dfMSE_b$RMSE,D)) #Misspecified 
  
  colnames(binary_t_table) <- c("Method","Abs bias","RMSE","Abs bias","RMSE","Abs bias","RMSE")
  binary_t_table
  
  write.csv(binary_t_table,file="data_results/binary_t_table.csv")
  
  binary_t_table2 <- data.frame(temp2_table_b,
                               round(dfBAL_b$Balance,D),round(dfTIME_b$Time,D)) #Misspecified 
  
  colnames(binary_t_table2) <- c("Method","Balance","TIME","Balance","TIME","Balance","TIME")
  binary_t_table2
  write.csv(binary_t_table2,file="data_results/binary_t_table2.csv")
  
  
  #########################################################################################################
  #######Continuous treatment - Table
  continuous_t_table <- data.frame(temp_table_c,
                               round(dfBIAS_c$`Abs-Bias`,D),round(dfMSE_c$RMSE,D)) #Misspecified 
  
  colnames(continuous_t_table) <- c("Method","Abs bias","RMSE","Abs bias","RMSE","Abs bias","RMSE")
  continuous_t_table
  
  write.csv(continuous_t_table,file="data_results/continuous_t_table.csv")
  continuous_t_table2 <- data.frame(temp2_table_c,
                                round(dfBAL_c$Balance,D),round(dfTIME_c$Time,D)) #Misspecified 
  
  colnames(continuous_t_table2) <- c("Method","Balance","TIME","Balance","TIME","Balance","TIME")
  continuous_t_table2
  
  write.csv(continuous_t_table2,file="data_results/continuous_t_table2.csv")
  
  

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
  
  D <- 3
  
  #########################################################################################################
  #########################################################################################################
  # Positivity
  #########################################################################################################
  #########################################################################################################
  
  load(paste("data_results/simu_errors_n",n,"positivity.Rdata"))
  
  #########################################################################################################
  #######SE
  dfSE_b <- dfSE[which(dfSE$Treatment=='binary' & dfSE$Type!="W-Bootstrap"),]
  dfSE_b$id <- sort(rep(1:(nrow(dfSE_b)/2),2))
  dfSE_b <- reshape(dfSE_b,timevar ="Positivity",direction="wide",idvar="id")
  
  dfSE_c <- dfSE[which(dfSE$Treatment=='continuous'  & dfSE$Type!="W-Bootstrap"),]
  dfSE_c$id <- sort(rep(1:(nrow(dfSE_c)/2),2))
  dfSE_c <- reshape(dfSE_c,timevar ="Positivity",direction="wide",idvar="id")
  
  
  #########################################################################################################
  #######COVERAGE
  dfCO_b <- dfCO[which(dfCO$Treatment=='binary' & dfCO$Type!="W-Bootstrap-Perc" & dfCO$Type!="W-Bootstrap-Norm"),]
  dfCO_b$id <- sort(rep(1:(nrow(dfCO_b)/2),2))
  dfCO_b <- reshape(dfCO_b,timevar ="Positivity",direction="wide",idvar="id")
  
  dfCO_c <- dfCO[which(dfCO$Treatment=='continuous' & dfCO$Type!="W-Bootstrap-Perc" & dfCO$Type!="W-Bootstrap-Norm"),]
  dfCO_c$id <- sort(rep(1:(nrow(dfCO_c)/2),2))
  dfCO_c <- reshape(dfCO_c,timevar ="Positivity",direction="wide",idvar="id")
  
  
  #########################################################################################################
  #######TIME
  dfTIME_b <- dfTIME[which(dfTIME$Treatment=='binary' & dfTIME$Type!="W-Bootstrap"),]
  dfTIME_b$id <- sort(rep(1:(nrow(dfTIME_b)/2),2))
  dfTIME_b <- reshape(dfTIME_b,timevar ="Positivity",direction="wide",idvar="id")
  
  
  dfTIME_c <- dfTIME[which(dfTIME$Treatment=='continuous' & dfTIME$Type!="W-Bootstrap"),]
  dfTIME_c$id <- sort(rep(1:(nrow(dfTIME_c)/2),2))
  dfTIME_c <- reshape(dfTIME_c,timevar ="Positivity",direction="wide",idvar="id")
  
  
  
  #########################################################################################################
  #######Binary treatment - Table
  temp_table_bn <- data.frame(c("Empirical","Bootstrap-Perc","Bootstrap-Norm","Robust","Naive"),
                             c(round(dfSE_b$`SE.1`,D)[1:2],round(dfSE_b$`SE.1`,D)[2],round(dfSE_b$`SE.1`,D)[3:4]),
                             c("-",round(dfCO_b$Coverage.1,D)),
                             c("-",round(dfTIME_b$`Time.1`,D)[1:2],round(dfTIME_b$`Time.1`,D)[1],round(dfTIME_b$`Time.1`,D)[3]), #Correct - no PPV
                             c(round(dfSE_b$`SE.2`,D)[1:2],round(dfSE_b$`SE.2`,D)[2],round(dfSE_b$`SE.2`,D)[3:4]),
                             c("-",round(dfCO_b$Coverage.2,D)),
                             c("-",round(dfTIME_b$`Time.2`,D)[1],round(dfTIME_b$`Time.2`,D)[1],round(dfTIME_b$`Time.2`,D)[2],round(dfTIME_b$`Time.2`,D)[3]) ) #Correct - PPV 
  
  
  
  #########################################################################################################
  #######Continuous treatment - Table
  temp_table_cn <- data.frame(c("Empirical","Bootstrap-Perc","Bootstrap-Norm","Robust","Naive"),
                             c(round(dfSE_c$`SE.1`,D)[1:2],round(dfSE_c$`SE.1`,D)[2],round(dfSE_c$`SE.1`,D)[3:4]),
                             c("-",round(dfCO_c$Coverage.1,D)),
                             c("-",round(dfTIME_c$`Time.1`,D)[1:2],round(dfTIME_c$`Time.1`,D)[1],round(dfTIME_c$`Time.1`,D)[3]), #Correct - no PPV
                             c(round(dfSE_c$`SE.2`,D)[1:2],round(dfSE_c$`SE.2`,D)[2],round(dfSE_c$`SE.2`,D)[3:4]),
                             c("-",round(dfCO_c$Coverage.2,D)),
                             c("-",round(dfTIME_c$`Time.2`,D)[1],round(dfTIME_c$`Time.2`,D)[1],round(dfTIME_c$`Time.2`,D)[1],round(dfTIME_c$`Time.2`,D)[3]) ) #Correct - PPV 
  
  
  
  
  
  
  #########################################################################################################
  #########################################################################################################
  # Misspecification
  #########################################################################################################
  #########################################################################################################
  
  load(paste("data_results/simu_errors_n",n,"misspecification.Rdata"))

  #########################################################################################################
  #######SE
  dfSE_b <- dfSE[which(dfSE$Treatment=='binary' & dfSE$Type!="W-Bootstrap"),]
  dfSE_c <- dfSE[which(dfSE$Treatment=='continuous' & dfSE$Type!="W-Bootstrap"),]
  
  #########################################################################################################
  #######SE
  dfCO_b <- dfCO[which(dfCO$Treatment=='binary' & dfCO$Type!="W-Bootstrap-Norm" & dfCO$Type!="W-Bootstrap-Perc"),]
  dfCO_c <- dfCO[which(dfCO$Treatment=='continuous' & dfCO$Type!="W-Bootstrap-Perc" & dfCO$Type!="W-Bootstrap-Norm"),]
  
  #########################################################################################################
  #######SE
  dfTIME_b <- dfTIME[which(dfTIME$Treatment=='binary' & dfTIME$Type!="W-Bootstrap"),]
  dfTIME_c <- dfTIME[which(dfTIME$Treatment=='continuous' & dfTIME$Type!="W-Bootstrap"),]
  
  
  
  #########################################################################################################
  #######Binary treatment - Table
  binary_t_table_row <- data.frame( temp_table_bn,
                                    c(round(dfSE_b$`SE`,D)[1:2],round(dfSE_b$`SE`,D)[2],round(dfSE_b$`SE`,D)[3:4]),
                                    c("-",round(dfCO_b$Coverage,D)),
                                    c("-",round(dfTIME_b$`Time`,D)[1],round(dfTIME_b$`Time`,D)[1],round(dfTIME_b$`Time`,D)[1],round(dfTIME_b$`Time`,D)[3])
                                    ) #Misspecified 
  
  colnames(binary_t_table_row) <- c("Method","SE","Coverage","Time","SE","Coverage","Time","SE","Coverage","Time")
  binary_t_table_row
  
  write.csv(binary_t_table_row,file="data_results/binary_t_table_row.csv")
  
  #########################################################################################################
  #######Continuous treatment - Table
  continuous_t_table_row <- data.frame( temp_table_cn,
                                        c(round(dfSE_c$`SE`,D)[1:2],round(dfSE_c$`SE`,D)[2],round(dfSE_c$`SE`,D)[3:4]),
                                        c("-",round(dfCO_c$Coverage,D)),
                                        c("-",round(dfTIME_c$`Time`,D)[1],round(dfTIME_c$`Time`,D)[1],round(dfTIME_c$`Time`,D)[1],round(dfTIME_c$`Time`,D)[3])
  ) #Misspecified 
  
  colnames(continuous_t_table_row) <- c("Method","SE","Coverage","Time","SE","Coverage","Time","SE","Coverage","Time")
  continuous_t_table_row
  
  write.csv(continuous_t_table_row,file="data_results/continuous_t_table_row.csv")
  
  

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

# dfBIAS_b <- dfBIAS[which(dfBIAS$Treatment=='binary') ,]
# dfBIAS_c <- dfBIAS[which(dfBIAS$Treatment=='continuous'),]

gBiasD_b <- ggplot(dfBIAS_b, aes(x=(Delta),y=(abs(Bias)), col=Method)) +
  geom_line(aes(linetype=Method),size=1) +
  geom_point(aes(linetype=Method),size=3) +
  xlab(expression(paste(delta))) +
  scale_x_continuous(breaks=c(0, 0.05, 0.1),labels=c("0.00001", "0.05", "0.1")) +
  coord_cartesian(ylim = c(0, 0.8)) +
  scale_color_manual(values=c("black")) +
  theme_bw() + theme(legend.position="none") +
  ylab("Abs(Bias)")


gBiasD_c <- ggplot(dfBIAS_c, aes(x=Delta,y=abs(Bias), col=Method)) +
  geom_line(aes(linetype=Method),size=1) +
  geom_point(aes(linetype=Method),size=3) +
  xlab(expression(paste(delta))) +
  scale_x_continuous(breaks=c(0, 0.05, 0.1),labels=c("0.00001", "0.05", "0.1")) +
  coord_cartesian(ylim = c(0, 0.8)) +
  scale_color_manual(values=c("black")) +
  theme_bw() + theme(legend.position="none") +
  ylab("Abs(Bias)")



#########################################################################################################
#######MSE
dfMSE_b <- dfMSE[which(dfMSE$Treatment=='binary' & dfBIAS$Delta<=0.1),]
dfMSE_c <- dfMSE[which(dfMSE$Treatment=='continuous' & dfBIAS$Delta<=0.1),]

# dfMSE_b <- dfMSE[which(dfMSE$Treatment=='binary'),]
# dfMSE_c <- dfMSE[which(dfMSE$Treatment=='continuous'),]


gMseD_b <- ggplot(dfMSE_b, aes(x=Delta,y=sqrt(MSE), col=Method)) +
  geom_line(aes(linetype=Method),size=1) +
  geom_point(aes(linetype=Method),size=3) +
  xlab(expression(paste(delta))) +
  scale_x_continuous(breaks=c(0, 0.05, 0.1),labels=c("0.00001", "0.05", "0.1")) +
  coord_cartesian(ylim = c(0, 0.8)) +
  scale_color_manual(values=c("black")) +
  theme_bw() + theme(legend.position="none") +
  ylab("RMSE")


gMseD_c <- ggplot(dfMSE_c, aes(x=Delta,y=sqrt(MSE), col=Method)) +
  geom_line(aes(linetype=Method),size=1) +
  geom_point(aes(linetype=Method),size=3) +
  xlab(expression(paste(delta))) +
  scale_x_continuous(breaks=c(0, 0.05, 0.1),labels=c("0.00001", "0.05", "0.1")) +
  coord_cartesian(ylim = c(0, 0.8)) +
  scale_color_manual(values=c("black")) +
  theme_bw() + theme(legend.position="none") +
  ylab("RMSE")



#########################################################################################################
#######BAL
dfBAL_b <- dfBAL[which(dfBAL$Treatment=='binary' & dfBIAS$Delta<=0.1),]
dfBAL_c <- dfBAL[which(dfBAL$Treatment=='continuous' & dfBIAS$Delta<=0.1),]

# dfBAL_b <- dfBAL[which(dfBAL$Treatment=='binary'),]
# dfBAL_c <- dfBAL[which(dfBAL$Treatment=='continuous'),]


gBalD_b <- ggplot(dfBAL_b, aes(x=Delta,y=Balance, col=Method)) +
  geom_line(aes(linetype=Method),size=1) +
  geom_point(aes(linetype=Method),size=3) +
  xlab(expression(paste(delta))) +
  scale_x_continuous(breaks=c(0, 0.05, 0.1),labels=c("0.00001", "0.05", "0.1")) +
  coord_cartesian(ylim = c(0, 1)) +
  scale_color_manual(values=c("black")) +
  theme_bw() + theme(legend.position="none") +
  ylab("Balance")


gBalD_c <- ggplot(dfBAL_c, aes(x=Delta,y=Balance, col=Method)) +
  geom_line(aes(linetype=Method),size=1) +
  geom_point(aes(linetype=Method),size=3) +
  xlab(expression(paste(delta))) +
  scale_x_continuous(breaks=c(0, 0.05, 0.1),labels=c("0.00001", "0.05", "0.1")) +
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
  scale_x_continuous(breaks=c(0, 0.05, 0.1),labels=c("0.00001", "0.05", "0.1")) +
  coord_cartesian(ylim = c(0, 0.05)) +
  scale_color_manual(values=c("black")) +
  theme_bw() + theme(legend.position="none") +
  ylab("Computational time (sec)")


gTimeD_c <- ggplot(dfTIME_c, aes(x=Delta,y=Time, col=Method)) +
  geom_line(aes(linetype=Method),size=1) +
  geom_point(aes(linetype=Method),size=3) +
  xlab(expression(paste(delta))) +
  scale_x_continuous(breaks=c(0, 0.05, 0.1),labels=c("0.00001", "0.05", "0.1")) +
  coord_cartesian(ylim = c(0, 0.05)) +
  scale_color_manual(values=c("black")) +
  theme_bw() + theme(legend.position="none") +
  ylab("Computational time (sec)")




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


#For Reviewer
grid.newpage()
pushViewport(viewport(layout = grid.layout(2, 2, heights = unit(c(0.5,5,
                                                                  0.5,5), "null"))))
grid.text("Binary treatment", vp = viewport(layout.pos.row = 1, layout.pos.col = 1))
print(gTimeD_b, vp = viewport(layout.pos.row = 2, layout.pos.col = 1))


grid.text("Continuous treatment", vp = viewport(layout.pos.row = 1, layout.pos.col = 2))
print(gTimeD_c, vp = viewport(layout.pos.row = 2, layout.pos.col = 2))

