########################################################################
########################################################################
#***********************************************************************
#
#   In this file, we provide code to evaluate the relationship between
#   covariate and treatment (Figure 16 of the Supplementary Material)
# 
#***********************************************************************
########################################################################
########################################################################

rm(list=ls())
set.seed(1)


library(gurobi)
library(ggplot2)
library(cobalt)
library(LaplacesDemon)
library(grid)

row <- function(confounders,intervention,delta=0.01){
  
  Xs <- scale(cbind(confounders))
  trs <- as.numeric(scale(intervention))
  if(is.null(dim(confounders)[1])!=T){
    n <- dim(confounders)[1]
  }
  else{
    n <- length(confounders)
  }#end else
  
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

SIZET <- 14
SIZETi <- 16
LENGTHLe <- 1
n <- 5000

###################################################################################################################
#Continuous

#***********************************************************************************************
#***********************************************************************************************
#  
# Linear Dependence  
#
#***********************************************************************************************
#***********************************************************************************************  

#***********************************************************************************************
#simulations setting
x         <- rnorm(n,0,1)
t         <- 1*x + rnorm(n,0,1)
ti        <- rnorm(n,0,1)
y         <- 1*t + 2*x + rnorm(n,0,1)

Tr        <- t
X         <-  x
Y         <-  y

data <- data.frame(X,Tr,Y)
colnames(data) <- c("X","Tr","Y")


confounders <- data.frame(X)
intervention <- t

#***********************************************************************************************
#solving QP
row_w <- row(confounders,intervention,delta=0.0001) 


#***********************************************************************************************
#plot 
weights_plot_correct_c <- row_w$x*n

Xs <- scale(X)
Trs <- scale(Tr)
data    <- data.frame(Xs,Trs,Y,weights_plot_correct_c) 


p1 <- ggplot(data, aes(x = Xs, y = Trs,size=weights_plot_correct_c,fill=weights_plot_correct_c)) + 
  geom_point(shape = 21) + 
  geom_smooth(aes(x = Xs, y = Trs, weight=weights_plot_correct_c, col="red"),se=F,method="glm", method.args = list(family = "gaussian")) +
  geom_smooth(aes(x = Xs, y = Trs),se=F) +
  scale_color_manual(values=c("black")) + 
  scale_fill_continuous(low = "grey85", high = "black") +
  theme_bw() + labs(size = "ROW", fill = "ROW") +
  ggtitle("Linear Dependence") + 
  labs(x = "Covariate", y = "Treatment") +
  theme_bw() + theme(legend.position="none", 
                     axis.title = element_text(size = SIZET),
                     axis.text = element_text(size = SIZET),
                     plot.title = element_text(hjust = 0.5, size = SIZETi),
                     legend.key.width = unit(LENGTHLe,"cm")) +
  guides(color = FALSE, size = FALSE)  +
  ylim(min(Trs),max(Trs)) + xlim(min(Xs),max(Xs)) + 
  theme(legend.position="bottom",legend.direction="horizontal") +
  theme(legend.position="bottom",legend.direction="horizontal")



#***********************************************************************************************
#check results
bal.tab(Tr~X,weights=row_w$x,data=data)
summary(lm(y~t,weights = row_w$x))



#***********************************************************************************************
#***********************************************************************************************
#  
# Non-linear dependence with correlation  
#
#***********************************************************************************************
#***********************************************************************************************  

rm(x,t,ti,y,Tr,X,Y,data,confounders,intervention,weights_plot_correct_c,Xs,Trs)

#***********************************************************************************************
#simulations setting
x         <- rnorm(n,0,1)
t         <- x + 1*x^2 + rnorm(n,0,1)
ti        <- rnorm(n,0,1)
y         <- 1*t + 2*x + rnorm(n,0,1)

Tr        <- t
X         <-  x
Y         <-  y

data      <- data.frame(X,Tr,Y)
colnames(data) <- c("X","Tr","Y")

confounders <- data.frame(X, X^2)
intervention <- t

#***********************************************************************************************
#solving QP
row_w <- row(confounders,intervention,delta=0.0001) 


#***********************************************************************************************
#plot 
weights_plot_correct_c <- row_w$x*n

Xs <- scale(X)
Trs <- scale(Tr)
data    <- data.frame(Xs,Trs,Y,weights_plot_correct_c) 


p2 <- ggplot(data, aes(x = Xs, y = Trs,size=weights_plot_correct_c,fill=weights_plot_correct_c)) + 
  geom_point(shape = 21) + 
  geom_smooth(aes(x = Xs, y = Trs, weight=weights_plot_correct_c, col="red"),se=F,method="glm", method.args = list(family = "gaussian")) +
  geom_smooth(aes(x = Xs, y = Trs),se=F) +
  scale_color_manual(values=c("black")) + 
  scale_fill_continuous(low = "grey85", high = "black") +
  theme_bw() + labs(size = "ROW", fill = "ROW") +
  ggtitle("Nonlinear dependence (quadratic)") + 
  labs(x = "Covariate", y = "Treatment") +
  theme_bw() + theme(legend.position="none", 
                     axis.title = element_text(size = SIZET),
                     axis.text = element_text(size = SIZET),
                     plot.title = element_text(hjust = 0.5, size = SIZETi),
                     legend.key.width = unit(LENGTHLe,"cm")) +
  guides(color = FALSE, size = FALSE)  +
  ylim(min(Trs),max(Trs)) + xlim(min(Xs),max(Xs)) + 
  theme(legend.position="bottom",legend.direction="horizontal") +
  theme(legend.position="bottom",legend.direction="horizontal")



#***********************************************************************************************
#check results
bal.tab(Tr~X,weights=row_w$x,data=data)
summary(lm(y~t,weights = row_w$x))



#***********************************************************************************************
#***********************************************************************************************
#  
# Non-linear dependence without correlation  
#
#***********************************************************************************************
#***********************************************************************************************  

rm(x,t,ti,y,Tr,X,Y,data,confounders,intervention,weights_plot_correct_c,Xs,Trs)

#***********************************************************************************************
#simulations setting

x                               <- runif(n,min=-.5,max=.5)  
p                               <- c(0.5,0.5)
mu                              <- c(-1, 1)
sigma                           <- c(1/3,1/3)
t                               <- rep(NA,n)
t[which( (x<(-1/6)) )] <- rnormm(length( which( (x<(-1/6))  )), p, mu, sigma)
t[which( (x>1/6))] <- rnormm(length( which( (x>1/6)) ), p, mu, sigma)
t[which( (x>=-1/6) & (x<=1/6))] <- rnorm(length( which( (x>=(-1/6)) & (x<=1/6)) ),0,1/3)
y         <- 1*t + 2*x + rnorm(n,0,1)
Gr        <- c( rep("Correlated",n) )

Tr        <- t
X         <-  x
Y         <-  y

data      <- data.frame(X,Tr,Y)
colnames(data) <- c("X","Tr","Y")


confounders <- data.frame(X)
intervention <- t

#***********************************************************************************************
#solving QP
row_w <- row(confounders,intervention,delta=0.0001) 


#***********************************************************************************************
#plot 
weights_plot_correct_c <- row_w$x*n

Xs <- scale(X)
Trs <- scale(Tr)
data    <- data.frame(Xs,Trs,Y,weights_plot_correct_c) 

p3 <- ggplot(data, aes(x = Xs, y = Trs,size=weights_plot_correct_c,fill=weights_plot_correct_c)) + 
  geom_point(shape = 21) + 
  geom_smooth(aes(x = Xs, y = Trs, weight=weights_plot_correct_c, col="red"),se=F,method="glm", method.args = list(family = "gaussian")) +
  geom_smooth(aes(x = Xs, y = Trs),se=F) +
  scale_color_manual(values=c("black")) + 
  scale_fill_continuous(low = "grey85", high = "black") +
  theme_bw() + labs(size = "ROW", fill = "ROW") +
  ggtitle("Nonlinear dependence w/o correlation  ") + 
  labs(x = "Covariate", y = "Treatment") +
  theme_bw() + theme(legend.position="none", 
                     axis.title = element_text(size = SIZET),
                     axis.text = element_text(size = SIZET),
                     plot.title = element_text(hjust = 0.5, size = SIZETi),
                     legend.key.width = unit(LENGTHLe,"cm")) +
  guides(color = FALSE, size = FALSE)  +
  ylim(min(Trs),max(Trs)) + xlim(min(Xs),max(Xs)) + 
  theme(legend.position="bottom",legend.direction="horizontal") +
  theme(legend.position="bottom",legend.direction="horizontal")


#***********************************************************************************************
#check results
bal.tab(Tr~X,weights=row_w$x,data=data)
summary(lm(y~t,weights = row_w$x))



#***********************************************************************************************
#***********************************************************************************************
#  
#  Independence
#
#***********************************************************************************************
#***********************************************************************************************  

rm(x,t,ti,y,Tr,X,Y,data,confounders,intervention,weights_plot_correct_c,Xs,Trs)

#***********************************************************************************************
#simulations setting
x         <- rnorm(n,0,1)
t         <- rnorm(n,0,1)
y         <- 1*t + 2*x + rnorm(n,0,1)

Tr        <- t
X         <-  x
Y         <-  y

data      <- data.frame(X,Tr,Y)
colnames(data) <- c("X","Tr","Y")


confounders <- data.frame(X)
intervention <- t

#***********************************************************************************************
#solving QP
row_w <- row(confounders,intervention,delta=0.0001) 


#***********************************************************************************************
#plot 
weights_plot_correct_c <- row_w$x*n

Xs <- scale(X)
Trs <- scale(Tr)
data    <- data.frame(Xs,Trs,Y,weights_plot_correct_c) 

p4 <- ggplot(data, aes(x = Xs, y = Trs,size=weights_plot_correct_c,fill=weights_plot_correct_c)) + 
  geom_point(shape = 21) + 
  geom_smooth(aes(x = Xs, y = Trs, weight=weights_plot_correct_c, col="red"),se=F,method="glm", method.args = list(family = "gaussian")) +
  geom_smooth(aes(x = Xs, y = Trs),se=F) +
  scale_color_manual(values=c("black")) + 
  scale_fill_continuous(low = "grey85", high = "black") +
  theme_bw() + labs(size = "ROW", fill = "ROW") +
  ggtitle("Independence") + 
  labs(x = "Covariate", y = "Treatment") +
  theme_bw() + theme(legend.position="none", 
                     axis.title = element_text(size = SIZET),
                     axis.text = element_text(size = SIZET),
                     plot.title = element_text(hjust = 0.5, size = SIZETi),
                     legend.key.width = unit(LENGTHLe,"cm")) +
  guides(color = FALSE, size = FALSE)  +
  ylim(min(Trs),max(Trs)) + xlim(min(Xs),max(Xs)) + 
  theme(legend.position="bottom",legend.direction="horizontal") +
  theme(legend.position="bottom",legend.direction="horizontal")



#***********************************************************************************************
#check results
bal.tab(Tr~X,weights=row_w$x,data=data)
summary(lm(y~t,weights = row_w$x))







#***********************************************************************************************
#***********************************************************************************************
#  
#  Sin
#
#***********************************************************************************************
#***********************************************************************************************  

rm(x,t,ti,y,Tr,X,Y,data,confounders,intervention,weights_plot_correct_c,Xs,Trs)

#***********************************************************************************************
#simulations setting
x         <- rnorm(n,0,4)
t         <- sin(x) + rnorm(n,0,0.1)
y         <- 1*t + 2*x + rnorm(n,0,1)


Tr        <- t
X         <-  x
Y         <-  y

data      <- data.frame(X,Tr,Y)
colnames(data) <- c("X","Tr","Y")

confounders <- data.frame(X)
intervention <- t

#***********************************************************************************************
#solving QP
row_w <- row(confounders,intervention,delta=0.0001) 


#***********************************************************************************************
#plot 
weights_plot_correct_c <- row_w$x*n

Xs <- scale(X)
Trs <- scale(Tr)
data    <- data.frame(Xs,Trs,Y,weights_plot_correct_c) 

p5 <- ggplot(data, aes(x = Xs, y = Trs,size=weights_plot_correct_c,fill=weights_plot_correct_c)) + 
  geom_point(shape = 21) + 
  geom_smooth(aes(x = Xs, y = Trs, weight=weights_plot_correct_c, col="red"),se=F,method="glm", method.args = list(family = "gaussian")) +
  geom_smooth(aes(x = Xs, y = Trs),se=F) +
  scale_color_manual(values=c("black")) + 
  scale_fill_continuous(low = "grey85", high = "black") +
  theme_bw() + labs(size = "ROW", fill = "ROW") +
  ggtitle("Sinusoidal dependence") + 
  labs(x = "Covariate", y = "Treatment") +
  theme_bw() + theme(legend.position="none", 
                     axis.title = element_text(size = SIZET),
                     axis.text = element_text(size = SIZET),
                     plot.title = element_text(hjust = 0.5, size = SIZETi),
                     legend.key.width = unit(LENGTHLe,"cm")) +
  guides(color = FALSE, size = FALSE)  +
  ylim(min(Trs),max(Trs)) + xlim(min(Xs),max(Xs)) + 
  theme(legend.position="bottom",legend.direction="horizontal") +
  theme(legend.position="bottom",legend.direction="horizontal")



#***********************************************************************************************
#check results
bal.tab(Tr~X,weights=row_w$x,data=data)
summary(lm(y~t,weights = row_w$x))

# heatscatter(X,Tr,colpal="matlablike2")




#***********************************************************************************************
#***********************************************************************************************
#  
#  Higher polynomiual
#
#***********************************************************************************************
#***********************************************************************************************  

rm(x,t,ti,y,Tr,X,Y,data,confounders,intervention,weights_plot_correct_c,Xs,Trs)

#***********************************************************************************************
#simulations setting
# x         <- rnorm(n,1,1)
# t         <- 1/x + rnorm(n,0,5)
x         <- runif(n,-10,10)
t         <- 0.5*(x+0.1)^3 + rnorm(n,0,1)
y         <- 1*t + 2*x + rnorm(n,0,1)

Tr        <- t
X         <-  x
Y         <-  y

data      <- data.frame(X,Tr,Y)
colnames(data) <- c("X","Tr","Y")    

confounders <- data.frame(X,X^2,X^3)
intervention <- t

#***********************************************************************************************
#solving QP
row_w <- row(confounders,intervention,delta=0.0001) 


#***********************************************************************************************
#plot 
weights_plot_correct_c <- row_w$x*n

Xs <- scale(X)
Trs <- scale(Tr)
data    <- data.frame(Xs,Trs,Y,weights_plot_correct_c) 

p6 <- ggplot(data, aes(x = Xs, y = Trs,size=weights_plot_correct_c,fill=weights_plot_correct_c)) + 
  geom_point(shape = 21) + 
  geom_smooth(aes(x = Xs, y = Trs, weight=weights_plot_correct_c, col="red"),se=F,method="glm", method.args = list(family = "gaussian")) +
  geom_smooth(aes(x = Xs, y = Trs),se=F) +
  scale_color_manual(values=c("black")) + 
  scale_fill_continuous(low = "grey85", high = "black") +
  theme_bw() + labs(size = "ROW", fill = "ROW") +
  ggtitle("Nonlinear dependence (cubic)") + 
  labs(x = "Covariate", y = "Treatment") +
  theme_bw() + theme(legend.position="none", 
                     axis.title = element_text(size = SIZET),
                     axis.text = element_text(size = SIZET),
                     plot.title = element_text(hjust = 0.5, size = SIZETi),
                     legend.key.width = unit(LENGTHLe,"cm")) +
  guides(color = FALSE, size = FALSE)  +
  ylim(min(Trs),max(Trs)) + xlim(min(Xs),max(Xs)) + 
  theme(legend.position="bottom",legend.direction="horizontal") +
  theme(legend.position="bottom",legend.direction="horizontal")



#***********************************************************************************************
#check results
bal.tab(Tr~X,weights=row_w$x,data=data)
summary(lm(y~t,weights = row_w$x))






#***********************************************************************************************
#***********************************************************************************************
#  
#  Right-skewed treatment
#
#***********************************************************************************************
#***********************************************************************************************  

rm(x,t,ti,y,Tr,X,Y,data,confounders,intervention,weights_plot_correct_c,Xs,Trs)

#***********************************************************************************************
#simulations setting
# x         <- rnorm(n,1,1)
# t         <- 1/x + rnorm(n,0,5)
x         <- rbeta(n,1,5)
t         <- 4*x + rlnorm(n,0,0.7)

data <- data.frame(t)
colnames(data) <- "Treatment"

pd7 <- ggplot(data, aes(Treatment)) +
  geom_density() + theme_bw()

y         <- 1*t + 2*x + rnorm(n,0,1)

Tr        <- t
X         <-  x
Y         <-  y

data      <- data.frame(X,Tr,Y)
colnames(data) <- c("X","Tr","Y")    

confounders <- data.frame(X)
intervention <- t

#***********************************************************************************************
#solving QP
row_w <- row(confounders,intervention,delta=0.0001) 


#***********************************************************************************************
#plot 
weights_plot_correct_c <- row_w$x*n

Xs <- scale(X)
Trs <- scale(Tr)
data    <- data.frame(Xs,Trs,Y,weights_plot_correct_c) 

p7 <- ggplot(data, aes(x = Xs, y = Trs,size=weights_plot_correct_c,fill=weights_plot_correct_c)) + 
  geom_point(shape = 21) + 
  geom_smooth(aes(x = Xs, y = Trs, weight=weights_plot_correct_c, col="red"),se=F,method="glm", method.args = list(family = "gaussian")) +
  geom_smooth(aes(x = Xs, y = Trs),se=F) +
  scale_color_manual(values=c("black")) + 
  scale_fill_continuous(low = "grey85", high = "black") +
  theme_bw() + labs(size = "ROW", fill = "ROW") +
  ggtitle("Right-skewed treatment") + 
  labs(x = "Covariate", y = "Treatment") +
  theme_bw() + theme(legend.position="none", 
                     axis.title = element_text(size = SIZET),
                     axis.text = element_text(size = SIZET),
                     plot.title = element_text(hjust = 0.5, size = SIZETi),
                     legend.key.width = unit(LENGTHLe,"cm")) +
  guides(color = FALSE, size = FALSE)  +
  ylim(min(Trs),max(Trs)) + xlim(min(Xs),max(Xs)) + 
  theme(legend.position="bottom",legend.direction="horizontal") +
  theme(legend.position="bottom",legend.direction="horizontal")



#***********************************************************************************************
#check results
bal.tab(Tr~X,weights=row_w$x,data=data)
summary(lm(y~t,weights = row_w$x))



#***********************************************************************************************
#***********************************************************************************************
#  
#  Left-skewed treatment
#
#***********************************************************************************************
#***********************************************************************************************  

rm(x,t,ti,y,Tr,X,Y,data,confounders,intervention,weights_plot_correct_c,Xs,Trs)

#***********************************************************************************************
#simulations setting
x         <- rbeta(n,5,1)
t         <- 4*x + 5*rbeta(n,5,1)

data <- data.frame(t)
colnames(data) <- "Treatment"

pd8 <- ggplot(data, aes(Treatment)) +
  geom_density() + theme_bw()

y         <- 1*t + 2*x + rnorm(n,0,1)

Tr        <- t
X         <-  x
Y         <-  y

data      <- data.frame(X,Tr,Y)
colnames(data) <- c("X","Tr","Y")    

confounders <- data.frame(X)
intervention <- t

#***********************************************************************************************
#solving QP
row_w <- row(confounders,intervention,delta=0.0001) 


#***********************************************************************************************
#plot 
weights_plot_correct_c <- row_w$x*n

Xs <- scale(X)
Trs <- scale(Tr)
data    <- data.frame(Xs,Trs,Y,weights_plot_correct_c) 

p8 <- ggplot(data, aes(x = Xs, y = Trs,size=weights_plot_correct_c,fill=weights_plot_correct_c)) + 
  geom_point(shape = 21) + 
  geom_smooth(aes(x = Xs, y = Trs, weight=weights_plot_correct_c, col="red"),se=F,method="glm", method.args = list(family = "gaussian")) +
  geom_smooth(aes(x = Xs, y = Trs),se=F) +
  scale_color_manual(values=c("black")) + 
  scale_fill_continuous(low = "grey85", high = "black") +
  theme_bw() + labs(size = "ROW", fill = "ROW") +
  ggtitle("Left-skewed treatment") + 
  labs(x = "Covariate", y = "Treatment") +
  theme_bw() + theme(legend.position="none", 
                     axis.title = element_text(size = SIZET),
                     axis.text = element_text(size = SIZET),
                     plot.title = element_text(hjust = 0.5, size = SIZETi),
                     legend.key.width = unit(LENGTHLe,"cm")) +
  guides(color = FALSE, size = FALSE)  +
  ylim(min(Trs),max(Trs)) + xlim(min(Xs),max(Xs)) + 
  theme(legend.position="bottom",legend.direction="horizontal") +
  theme(legend.position="bottom",legend.direction="horizontal")



#***********************************************************************************************
#check results
bal.tab(Tr~X,weights=row_w$x,data=data)
summary(lm(y~t,weights = row_w$x))



#########################################################################################################
#######


grid.newpage()
pushViewport(viewport(layout = grid.layout(12, 2, heights = unit(c(0, 0.25, 5, 
                                                                   0, 0.25,5, 
                                                                   0, 0.25,5, 
                                                                   0, 0.25,5, 
                                                                   0), "null"))))  


print(p1, vp = viewport(layout.pos.row = 3, layout.pos.col = 1))
print(p2, vp = viewport(layout.pos.row = 3, layout.pos.col = 2))

print(p6, vp = viewport(layout.pos.row = 6, layout.pos.col = 1))
print(p3, vp = viewport(layout.pos.row = 6, layout.pos.col = 2))

print(p5, vp = viewport(layout.pos.row = 9, layout.pos.col = 1))
print(p4, vp = viewport(layout.pos.row = 9, layout.pos.col = 2))

print(p7, vp = viewport(layout.pos.row = 12, layout.pos.col = 1))
print(p8, vp = viewport(layout.pos.row = 12, layout.pos.col = 2))











