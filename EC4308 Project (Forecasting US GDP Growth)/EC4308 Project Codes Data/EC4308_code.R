#remove all objects from memory
rm(list = ls())

#install.packages("ISLR")
#install.packages("glmnet")
#install.packages("hdm")
#install.packages("readxl")
#install.packages("pls")
#install.packages("githubinstall")
#install.packages("hdm")
#install.packages("sandwich")
#install.packages("randomForest")
#githubinstall("HDeconometrics")
#install.packages("dplyr") 
#install.packages("ggplot2") 
#install.packages("patchwork") 
library(readr)
library(dplyr)
library(ISLR)
library(glmnet)
library(hdm)
#library(readxl)
library(HDeconometrics)
library(githubinstall)
library(sandwich) #library to estimate variance for DM test regression using NeweyWest()
library(randomForest)
library(ggplot2)
library(patchwork)

df <- read.csv("final_cleaned.csv")[c(3:732),]
df2 <- df[-1] %>% relocate(INDPRO)
rownames(df2) <- df$sasdate

###############################################################
##Parameter Defining: Train/Test/Validation. Function Defining
###############################################################

#create function for MSE
MSE <- function(pred, truth){ #start and end body of the function by { } - same as a loop 
  return(mean((truth - pred)^2)) #end function with a return(output) statement. Here we can go straight to return because the object of interest is a simple function of inputs
}

# auxiliary function to compute root MSE:
RMSE <- function(pred, truth){ #start and end body of the function by { } - same as a loop 
  return(sqrt(mean((truth - pred)^2)))
}

### preliminary data manipulation 
# remove 'Week' variable, Y is the matrix of all variables (dependent variable and predictors)
Y <- as.matrix(df2[,])
#get the y variable - INDPRO
yy=Y[,1] 
#number of out-of-sample observations (test window)
nprev=243
#auxiliary:get the out-of-sample true values (last 180 obs of INDPRO)
oosy=tail(yy,nprev) 


set.seed(123)

x = Y[,-1]
comps.all=princomp(scale(x,scale=TRUE)) 
summary(comps.all)
eigenvalues=comps.all$sdev^2
plot(eigenvalues/sum(eigenvalues),type='b')
varianceshares = eigenvalues/sum(eigenvalues)
sum(varianceshares[1:10])

################################################################################
### function calls to run ar
################################################################################
source("func-ar.R")

bar1c=ar.rolling.window(Y,nprev,1,1,type="fixed") #1-step AR forecast
bar3c=ar.rolling.window(Y,nprev,1,3,type="fixed") #4-step AR forecast
bar6c=ar.rolling.window(Y,nprev,1,6,type="fixed") #8-step AR forecast
bar12c=ar.rolling.window(Y,nprev,1,12,type="fixed") #12-step AR forecast

### AR forecasts RMSE:
ar.rmse1=bar1c$errors[1]
ar.rmse3=bar3c$errors[1]
ar.rmse6=bar6c$errors[1]
ar.rmse12=bar12c$errors[1]

#Plot benchmark coefficients
arcoef.ts=ts(bar1c$coef, start=c(1999,9), end=c(2019,11), freq=12)
colnames(arcoef.ts)=c("Constant", "Phi1", "Phi2", "Phi3", "Phi4")
plot.ts(arcoef.ts, main="AR regression coefficients", cex.axis=1.4)

#AR forecasts RMSE:

ar.rmse1=bar1c$errors[1]
ar.rmse3=bar3c$errors[1]
ar.rmse6=bar6c$errors[1]
ar.rmse12=bar12c$errors[1]

#Create ts objects out of 1-step and 12-step benchmark forecasts
bench1.ts=ts(cbind(bar1c$pred,oosy), start=c(1999,9), end=c(2019,12), freq=12)
colnames(bench1.ts)=c("AR(4)","True Value")

bench12.ts=ts(cbind(bar12c$pred,oosy), start=c(1999,9), end=c(2019,12), freq=12)
colnames(bench1.ts)=c("AR(4)","True Value")

#Plot 1-step forecasts:

plot.ts(bench1.ts[,1], main="1-step Benchmark forecasts", cex.axis=1.5, lwd=1.8, col="red", ylab="INDPRO Growth", ylim = range(bench1.ts))
points(bench1.ts[,2], type="l", col="black",lwd=2)
legend("bottomleft",legend=c("AR(4)","INDPRO Growth"), col=c("red","black"), lty = 1)

#Plot 12-step forecasts:

plot.ts(bench12.ts[,1], main="12-step Benchmark forecasts", cex.axis=1.5, lwd=1.8, col="red", ylab="INDPRO Growth", ylim = range(bench1.ts))
points(bench12.ts[,2], type="l", col="black",lwd=2)
legend("bottomleft",legend=c("AR(4)","INDPRO Growth"), col=c("red", "black"), lty = 1)



############################################################
##LASSO + Post LASSO
############################################################

source("func-lasso-noPCADum.R")

alpha=1 #set alpha=1 for LASSO

#lasso bic: run forecasts and get rmse
lasso1c_bic=lasso.rolling.window(Y,nprev,1,1,alpha,IC="bic")
lasso3c_bic=lasso.rolling.window(Y,nprev,1,3,alpha,IC="bic")
lasso6c_bic=lasso.rolling.window(Y,nprev,1,6,alpha,IC="bic")
lasso12c_bic=lasso.rolling.window(Y,nprev,1,12,alpha,IC="bic")
lasso.rmse1_bic=lasso1c_bic$errors[1]
lasso.rmse3_bic=lasso3c_bic$errors[1]
lasso.rmse6_bic=lasso6c_bic$errors[1]
lasso.rmse12_bic=lasso12c_bic$errors[1]

#lasso aic: run forecasts and get rmse
lasso1c_aic=lasso.rolling.window(Y,nprev,1,1,alpha,IC="aic")
lasso3c_aic=lasso.rolling.window(Y,nprev,1,3,alpha,IC="aic")
lasso6c_aic=lasso.rolling.window(Y,nprev,1,6,alpha,IC="aic")
lasso12c_aic=lasso.rolling.window(Y,nprev,1,12,alpha,IC="aic")
lasso.rmse1_aic=lasso1c_aic$errors[1]
lasso.rmse3_aic=lasso3c_aic$errors[1]
lasso.rmse6_aic=lasso6c_aic$errors[1]
lasso.rmse12_aic=lasso12c_aic$errors[1]

#lasso aicc: run forecasts and get rmse
lasso1c_aicc=lasso.rolling.window(Y,nprev,1,1,alpha,IC="aicc")
lasso3c_aicc=lasso.rolling.window(Y,nprev,1,3,alpha,IC="aicc")
lasso6c_aicc=lasso.rolling.window(Y,nprev,1,6,alpha,IC="aicc")
lasso12c_aicc=lasso.rolling.window(Y,nprev,1,12,alpha,IC="aicc")
lasso.rmse1_aicc=lasso1c_aicc$errors[1]
lasso.rmse3_aicc=lasso3c_aicc$errors[1]
lasso.rmse6_aicc=lasso6c_aicc$errors[1]
lasso.rmse12_aicc=lasso12c_aicc$errors[1]

### post lasso using lasso selected on bic
pols.lasso1c_bic=pols.rolling.window(Y,nprev,1,1,lasso1c_bic$coef)
pols.lasso3c_bic=pols.rolling.window(Y,nprev,1,3,lasso3c_bic$coef)
pols.lasso6c_bic=pols.rolling.window(Y,nprev,1,6,lasso6c_bic$coef)
pols.lasso12c_bic=pols.rolling.window(Y,nprev,1,12,lasso12c_bic$coef)
plasso.rmse1_bic=pols.lasso1c_bic$errors[1]
plasso.rmse3_bic=pols.lasso3c_bic$errors[1]
plasso.rmse6_bic=pols.lasso6c_bic$errors[1]
plasso.rmse12_bic=pols.lasso12c_bic$errors[1]



############################################################
##Ridge
############################################################


############################################################
##Elastic Net
############################################################
alpha=0.5 #set alpha to 0.5 for elastic net

# elastic net bic: run forecasts and get rmse 
elnet1c_bic=lasso.rolling.window(Y,nprev,1,1,0.5,IC="bic")
elnet3c_bic=lasso.rolling.window(Y,nprev,1,3,0.5,IC="bic")
elnet6c_bic=lasso.rolling.window(Y,nprev,1,6,0.5,IC="bic")
elnet12c_bic=lasso.rolling.window(Y,nprev,1,12,0.5,IC="bic")
elnet.rmse1_bic=elnet1c_bic$errors[1]
elnet.rmse3_bic=elnet3c_bic$errors[1]
elnet.rmse6_bic=elnet6c_bic$errors[1]
elnet.rmse12_bic=elnet12c_bic$errors[1]


#######################################################
##Sparsity analysis over time plots for LASSO and ElNet
#######################################################

#Get nonzero coefficient numbers for different horizons (LASSO(BIC))
c1c=rowSums(lasso1c_bic$coef != 0)
c3c=rowSums(lasso3c_bic$coef != 0)
c6c=rowSums(lasso6c_bic$coef != 0)
c12c=rowSums(lasso12c_bic$coef != 0)

#Create a ts object for the plot
lcoef.ts=ts(cbind(c1c,c3c,c6c,c12c), start=c(1999,9), end=c(2019,12), freq=12)
colnames(lcoef.ts)=c("H1","H3","H6","H12")
#Plot numbers of nonzero coefficients across the test window
plot.ts(lcoef.ts, main="Sparsity Analysis for LASSO",cex.axis=1.5)

#Get nonzero coefficient numbers for different horizons (Elastic Net)
ce1c=rowSums(elnet1c_bic$coef != 0)
ce3c=rowSums(elnet3c_bic$coef != 0)
ce6c=rowSums(elnet6c_bic$coef != 0)
ce12c=rowSums(elnet12c_bic$coef != 0)

#Create a respective ts object for the plot
elcoef.ts=ts(cbind(ce1c,ce3c,ce6c,ce12c), start=c(1999,9), end=c(2019,12), freq=12)
colnames(elcoef.ts)=c("H1","H3","H6","H12")
#Plot numbers of nonzero coefficients across the test window:
plot.ts(elcoef.ts, main="Sparsity Analysis for Elastic Net",cex.axis=1.5)

################################################
##Plot ML forecasts for 1 and 12-steps (for visualization purpose only) AR LASSO
##############################################

#This is only for AR & LASSO, 1Step & 12Step
#Create the time series object collecting 1-step best=performing ML forecasts
ml1.ts=ts(cbind(oosy,bar1c$pred,lasso1c_bic$pred),  start=c(1999,9), end=c(2019,12), freq=12)
colnames(ml1.ts)=c("True Value","AR(4)","LASSO")

#Create the time series object collecting 1-step best=performing ML forecasts
ml12.ts=ts(cbind(oosy,bar12c$pred,lasso12c_bic$pred),  start=c(1999,9), end=c(2019,12), freq=12)
colnames(ml1.ts)=c("True Value","AR(4)","LASSO")

#Plot the graph for 1-step forecasts
plot.ts(ml1.ts[,1], main="1-step ML forecasts", cex.axis=1.5, lwd=1.8, ylab="INDPRO")
points(ml1.ts[,2], type="l", col="blue",lwd=1.8)
points(ml1.ts[,3], type="l", col="red",lwd=1.8)
legend("bottomleft", c("INDPRO Growth","AR(4)","LASSO"), lty=c(1,1,1) ,col=c("black","blue","red"))


#Plot the graph for 12-step forecasts
plot.ts(ml12.ts[,1], main="12-step ML forecasts", cex.axis=1.5, lwd=1.8, ylab="INDPRO")
points(ml12.ts[,2], type="l", col="blue",lwd=1.8)
points(ml12.ts[,3], type="l", col="red",lwd=1.8)
legend("bottomleft", c("INDPRO Growth","AR(4)","LASSO"), lty=c(1,1,1) ,col=c("black","blue","red"))


################################################
##Plot ML forecasts for 1 and 12-steps (for visualization purpose only) PLASSO ELNET
##############################################

#This is only for POST-LASSO & ElNet, 1Step & 12Step
#Create the time series object collecting 1-step best=performing ML forecasts
ml1.ts=ts(cbind(oosy,bar1c$pred,pols.lasso1c_bic$pred,elnet1c_bic$pred),  start=c(1999,9), end=c(2019,12), freq=12)
colnames(ml1.ts)=c("True Value","AR(4)","PLASSO","Elastic Net")

#Create the time series object collecting 12-step best=performing ML forecasts
ml12.ts=ts(cbind(oosy,bar1c$pred,pols.lasso12c_bic$pred,elnet12c_bic$pred),  start=c(1999,9), end=c(2019,12), freq=12)
colnames(ml1.ts)=c("True Value","AR(4)","PLASSO","Elastic Net")

#Plot the graph for 1-step forecasts
plot.ts(ml1.ts[,1], main="1-step ML forecasts", cex.axis=1.5, lwd=1.8, ylab="INDPRO")
points(ml1.ts[,2], type="l", col="blue",lwd=1.8)
points(ml1.ts[,3], type="l", col="red",lwd=1.8)
points(ml1.ts[,4], type="l", col="green",lwd=1.8)
legend("bottomleft", c("INDPRO Growth","AR(4)","PLASSO","Elastic Net"), lty=c(1,1,1) ,col=c("black","blue","red"))


#Plot the graph for 12-step forecasts
plot.ts(ml12.ts[,1], main="12-step ML forecasts", cex.axis=1.5, lwd=1.8, ylab="INDPRO")
points(ml12.ts[,2], type="l", col="blue",lwd=1.8)
points(ml12.ts[,3], type="l", col="red",lwd=1.8)
points(ml1.ts[,4], type="l", col="green",lwd=1.8)
legend("bottomleft", c("INDPRO Growth","AR(4)","PLASSO","Elastic Net"), lty=c(1,1,1) ,col=c("black","blue","red"))

################################################################################
### RANDOM FOREST. CALLS THE randomForestCode.R
################################################################################
source("randomForestCode.R")



################################################################################
### simple average with RF
################################################################################
avgpred1=(lasso1c_bic$pred+pols.lasso1c_bic$pred+elnet1c_bic$pred+pred_rf_one+bar1c$pred)/5
avgrmse1=sqrt(mean((oosy-avgpred1)^2))

avgpred3=(lasso3c_bic$pred+pols.lasso3c_bic$pred+elnet3c_bic$pred+pred_rf_three+bar3c$pred)/5
avgrmse3=sqrt(mean((oosy-avgpred3)^2))

avgpred6=(lasso6c_bic$pred+pols.lasso6c_bic$pred+elnet6c_bic$pred+pred_rf_six+bar6c$pred)/5
avgrmse6=sqrt(mean((oosy-avgpred6)^2))

avgpred12=(lasso12c_bic$pred+pols.lasso12c_bic$pred+elnet12c_bic$pred+pred_rf_twelve+bar12c$pred)/5
avgrmse12=sqrt(mean((oosy-avgpred12)^2))

################################################################################
### simple average without RF
################################################################################
avgpred_1=(lasso1c_bic$pred+pols.lasso1c_bic$pred+elnet1c_bic$pred+bar1c$pred)/4
avgrmse_1=sqrt(mean((oosy-avgpred_1)^2))

avgpred_3=(lasso3c_bic$pred+pols.lasso3c_bic$pred+elnet3c_bic$pred+bar3c$pred)/4
avgrmse_3=sqrt(mean((oosy-avgpred_3)^2))

avgpred_6=(lasso6c_bic$pred+pols.lasso6c_bic$pred+elnet6c_bic$pred+bar6c$pred)/4
avgrmse_6=sqrt(mean((oosy-avgpred_6)^2))

avgpred_12=(lasso12c_bic$pred+pols.lasso12c_bic$pred+elnet12c_bic$pred+bar12c$pred)/4
avgrmse_12=sqrt(mean((oosy-avgpred_12)^2))

############################################################
##Diebold Mariano Test on RMSE
############################################################
################################################################################
### DM test: AR(4) vs lasso 
################################################################################

# AR(4) squared loss for different horizons
loss_ar1c = (oosy-bar1c$pred)^2
loss_ar3c = (oosy-bar3c$pred)^2
loss_ar6c = (oosy-bar6c$pred)^2
loss_ar12c = (oosy-bar12c$pred)^2

# lasso : squared loss for different horizons
loss_lasso1c = (oosy-lasso1c_bic$pred)^2
loss_lasso3c = (oosy-lasso3c_bic$pred)^2
loss_lasso6c = (oosy-lasso6c_bic$pred)^2
loss_lasso12c = (oosy-lasso12c_bic$pred)^2

# compute loss differentials (d_t)
d_ar_lasso1 = loss_ar1c - loss_lasso1c
d_ar_lasso3 = loss_ar3c - loss_lasso3c
d_ar_lasso6 = loss_ar6c - loss_lasso6c
d_ar_lasso12 = loss_ar12c - loss_lasso12c

# check if stationarity assumption holds 
dt_ar_lasso.ts=ts(cbind(d_ar_lasso1,d_ar_lasso3,d_ar_lasso6,d_ar_lasso12), start=c(1999,9), end=c(2019,12), freq=12)
colnames(dt_ar_lasso.ts)=c("1-step dt","3-step dt","6-step dt","12-step dt")
plot.ts(dt_ar_lasso.ts, main="Loss differential",cex.axis=1.8)

#Regress d_t (AR-lasso) for 1-step forecasts on a constant - get estimate of mean(d_t)
dm_ar_lasso1=lm(d_ar_lasso1~1)
acf(dm_ar_lasso1$residuals)
dm_ar_lasso1$coefficients/sqrt(NeweyWest(dm_ar_lasso1,lag=6))

#3-step forecast test
dm_ar_lasso3=lm(d_ar_lasso3~1)
acf(dm_ar_lasso3$residuals)
dm_ar_lasso3$coefficients/sqrt(NeweyWest(dm_ar_lasso3,lag=6))

#6-step forecast test
dm_ar_lasso6=lm(d_ar_lasso6~1)
acf(dm_ar_lasso6$residuals)
dm_ar_lasso6$coefficients/sqrt(NeweyWest(dm_ar_lasso6,lag=6))

#12-step forecast test
dm_ar_lasso12=lm(d_ar_lasso12~1)
acf(dm_ar_lasso12$residuals)
dm_ar_lasso12$coefficients/sqrt(NeweyWest(dm_ar_lasso12,lag=6))

################################################################################
### DM test: AR(4) vs plasso 
################################################################################

# AR(4) squared loss for different horizons
loss_ar1c = (oosy-bar1c$pred)^2
loss_ar3c = (oosy-bar3c$pred)^2
loss_ar6c = (oosy-bar6c$pred)^2
loss_ar12c = (oosy-bar12c$pred)^2

# plasso : squared loss for different horizons
loss_plasso1c = (oosy-pols.lasso1c_bic$pred)^2
loss_plasso3c = (oosy-pols.lasso3c_bic$pred)^2
loss_plasso6c = (oosy-pols.lasso6c_bic$pred)^2
loss_plasso12c = (oosy-pols.lasso12c_bic$pred)^2

# compute loss differentials (d_t)
d_ar_plasso1 = loss_ar1c - loss_plasso1c
d_ar_plasso3 = loss_ar3c - loss_plasso3c
d_ar_plasso6 = loss_ar6c - loss_plasso6c
d_ar_plasso12 = loss_ar12c - loss_plasso12c

# check if stationarity assumption holds 
dt_ar_plasso.ts=ts(cbind(d_ar_plasso1,d_ar_plasso3,d_ar_plasso6,d_ar_plasso12), start=c(1999,9), end=c(2019,12), freq=12)
colnames(dt_ar_plasso.ts)=c("1-step dt","3-step dt","6-step dt","12-step dt")
plot.ts(dt_ar_plasso.ts, main="Loss differential",cex.axis=1.8)

#Regress d_t (AR-plasso) for 1-step forecasts on a constant - get estimate of mean(d_t)
dm_ar_plasso1=lm(d_ar_plasso1~1)
acf(dm_ar_plasso1$residuals)
dm_ar_plasso1$coefficients/sqrt(NeweyWest(dm_ar_plasso1,lag=6))

#3-step forecast test
dm_ar_plasso3=lm(d_ar_plasso3~1)
acf(dm_ar_plasso3$residuals)
dm_ar_plasso3$coefficients/sqrt(NeweyWest(dm_ar_plasso3,lag=6))

#6-step forecast test
dm_ar_plasso6=lm(d_ar_plasso6~1)
acf(dm_ar_plasso6$residuals)
dm_ar_plasso6$coefficients/sqrt(NeweyWest(dm_ar_plasso6,lag=6))

#12-step forecast test
dm_ar_plasso12=lm(d_ar_plasso12~1)
acf(dm_ar_plasso12$residuals)
dm_ar_plasso12$coefficients/sqrt(NeweyWest(dm_ar_plasso12,lag=6))


################################################################################
### DM test: AR(4) vs ElNet 
################################################################################

# AR(4) squared loss for different horizons
loss_ar1c = (oosy-bar1c$pred)^2
loss_ar3c = (oosy-bar3c$pred)^2
loss_ar6c = (oosy-bar6c$pred)^2
loss_ar12c = (oosy-bar12c$pred)^2

# elnet : squared loss for different horizons
loss_elnet1c = (oosy-elnet1c_bic$pred)^2
loss_elnet3c = (oosy-elnet3c_bic$pred)^2
loss_elnet6c = (oosy-elnet6c_bic$pred)^2
loss_elnet12c = (oosy-elnet12c_bic$pred)^2

# compute loss differentials (d_t)
d_ar_elnet1 = loss_ar1c - loss_elnet1c
d_ar_elnet3 = loss_ar3c - loss_elnet3c
d_ar_elnet6 = loss_ar6c - loss_elnet6c
d_ar_elnet12 = loss_ar12c - loss_elnet12c

# check if stationarity assumption holds 
dt_ar_elnet.ts=ts(cbind(d_ar_elnet1,d_ar_elnet3,d_ar_elnet6,d_ar_elnet12), start=c(1999,9), end=c(2019,12), freq=12)
colnames(dt_ar_elnet.ts)=c("H1 dt","H3 dt","H6 dt","H12 dt")
plot.ts(dt_ar_elnet.ts, main="Loss differential",cex.axis=1.8)

#Regress d_t (AR-elnet) for 1-step forecasts on a constant - get estimate of mean(d_t)
dm_ar_elnet1=lm(d_ar_elnet1~1)
acf(dm_ar_elnet1$residuals)
dm_ar_elnet1$coefficients/sqrt(NeweyWest(dm_ar_elnet1,lag=6))

#3-step forecast test
dm_ar_elnet3=lm(d_ar_elnet3~1)
acf(dm_ar_elnet3$residuals)
dm_ar_elnet3$coefficients/sqrt(NeweyWest(dm_ar_elnet3,lag=6))

#6-step forecast test
dm_ar_elnet6=lm(d_ar_elnet6~1)
acf(dm_ar_elnet6$residuals)
dm_ar_elnet6$coefficients/sqrt(NeweyWest(dm_ar_elnet6,lag=6))

#12-step forecast test
dm_ar_elnet12=lm(d_ar_elnet12~1)
acf(dm_ar_elnet12$residuals)
dm_ar_elnet12$coefficients/sqrt(NeweyWest(dm_ar_elnet12,lag=6))


################################################################################
### DM test: Elnet vs Plasso 
################################################################################

# elnet : squared loss for different horizons
loss_elnet1c = (oosy-elnet1c_bic$pred)^2
loss_elnet3c = (oosy-elnet3c_bic$pred)^2
loss_elnet6c = (oosy-elnet6c_bic$pred)^2
loss_elnet12c = (oosy-elnet12c_bic$pred)^2

# plasso : squared loss for different horizons
loss_plasso1c = (oosy-pols.lasso1c_bic$pred)^2
loss_plasso3c = (oosy-pols.lasso3c_bic$pred)^2
loss_plasso6c = (oosy-pols.lasso6c_bic$pred)^2
loss_plasso12c = (oosy-pols.lasso12c_bic$pred)^2

# compute loss differentials (d_t)
d_elnet_plasso1 = loss_elnet1c - loss_plasso1c
d_elnet_plasso3 = loss_elnet3c - loss_plasso3c
d_elnet_plasso6 = loss_elnet6c - loss_plasso6c
d_elnet_plasso12 = loss_elnet12c - loss_plasso12c

# check if stationarity assumption holds 
dt_elnet_plasso.ts=ts(cbind(d_elnet_plasso1,d_elnet_plasso3,d_elnet_plasso6,d_elnet_plasso12), start=c(1999,9), end=c(2019,12), freq=12)
colnames(dt_elnet_plasso.ts)=c("1-step dt","3-step dt","6-step dt","12-step dt")
plot.ts(dt_elnet_plasso.ts, main="Loss differential",cex.axis=1.8)

#Regress d_t (elnet-plasso) for 1-step forecasts on a constant - get estimate of mean(d_t)
dm_elnet_plasso1=lm(d_elnet_plasso1~1)
acf(dm_elnet_plasso1$residuals)
dm_elnet_plasso1$coefficients/sqrt(NeweyWest(dm_elnet_plasso1,lag=6))

#3-step forecast test
dm_elnet_plasso3=lm(d_elnet_plasso3~1)
acf(dm_elnet_plasso3$residuals)
dm_elnet_plasso3$coefficients/sqrt(NeweyWest(dm_elnet_plasso3,lag=6))

#6-step forecast test
dm_elnet_plasso6=lm(d_elnet_plasso6~1)
acf(dm_elnet_plasso6$residuals)
dm_elnet_plasso6$coefficients/sqrt(NeweyWest(dm_elnet_plasso6,lag=6))

#12-step forecast test
dm_elnet_plasso12=lm(d_elnet_plasso12~1)
acf(dm_elnet_plasso12$residuals)
dm_elnet_plasso12$coefficients/sqrt(NeweyWest(dm_elnet_plasso12,lag=6))


################################################################################
### DM test: AR(4) vs RF
################################################################################

# AR(4) squared loss for different horizons
loss_ar1c = (oosy-bar1c$pred)^2
loss_ar3c = (oosy-bar3c$pred)^2
loss_ar6c = (oosy-bar6c$pred)^2
loss_ar12c = (oosy-bar12c$pred)^2

# RF : squared loss for different horizons
loss_rf1c = (oosy-pred_rf_one)^2
loss_rf3c = (oosy-pred_rf_three)^2
loss_rf6c = (oosy-pred_rf_six)^2
loss_rf12c = (oosy-pred_rf_twelve)^2

# compute loss differentials (d_t)
d_ar_rf1 = loss_ar1c - loss_rf1c
d_ar_rf3 = loss_ar3c - loss_rf3c
d_ar_rf6 = loss_ar6c - loss_rf6c
d_ar_rf12 = loss_ar12c - loss_rf12c

# check if stationarity assumption holds 
dt_ar_rf.ts=ts(cbind(d_ar_rf1,d_ar_rf3,d_ar_rf6,d_ar_rf12), start=c(1999,9), end=c(2019,12), freq=12)
colnames(dt_ar_rf.ts)=c("1-step dt","3-step dt","6-step dt","12-step dt")
plot.ts(dt_ar_rf.ts, main="Loss differential",cex.axis=1.8)

#Regress d_t (AR-RF) for 1-step forecasts on a constant - get estimate of mean(d_t)
dm_ar_rf1=lm(d_ar_rf1~1)
acf(dm_ar_rf1$residuals)
dm_ar_rf1$coefficients/sqrt(NeweyWest(dm_ar_rf1,lag=6))

#3-step forecast test
dm_ar_rf3=lm(d_ar_rf3~1)
acf(dm_ar_rf3$residuals)
dm_ar_rf3$coefficients/sqrt(NeweyWest(dm_ar_rf3,lag=6))

#6-step forecast test
dm_ar_rf6=lm(d_ar_rf6~1)
acf(dm_ar_rf6$residuals)
dm_ar_rf6$coefficients/sqrt(NeweyWest(dm_ar_rf6,lag=6))

#12-step forecast test
dm_ar_rf12=lm(d_ar_rf12~1)
acf(dm_ar_rf12$residuals)
dm_ar_rf12$coefficients/sqrt(NeweyWest(dm_ar_rf12,lag=6))

################################################################################
### DM test: AR(4) vs Bagging
################################################################################

# AR(4) squared loss for different horizons
loss_ar1c = (oosy-bar1c$pred)^2
loss_ar3c = (oosy-bar3c$pred)^2
loss_ar6c = (oosy-bar6c$pred)^2
loss_ar12c = (oosy-bar12c$pred)^2

# bg : squared loss for different horizons
loss_bg1c = (oosy-pred_bagging_one)^2
loss_bg3c = (oosy-pred_bagging_three)^2
loss_bg6c = (oosy-pred_bagging_six)^2
loss_bg12c = (oosy-pred_bagging_twelve)^2

# compute loss differentials (d_t)
d_ar_bg1 = loss_ar1c - loss_bg1c
d_ar_bg3 = loss_ar3c - loss_bg3c
d_ar_bg6 = loss_ar6c - loss_bg6c
d_ar_bg12 = loss_ar12c - loss_bg12c

# check if stationarity assumption holds 
dt_ar_bg.ts=ts(cbind(d_ar_bg1,d_ar_bg3,d_ar_bg6,d_ar_bg12), start=c(1999,9), end=c(2019,12), freq=12)
colnames(dt_ar_bg.ts)=c("1-step dt","3-step dt","6-step dt","12-step dt")
plot.ts(dt_ar_bg.ts, main="Loss differential",cex.axis=1.8)

#Regress d_t (AR-bg) for 1-step forecasts on a constant - get estimate of mean(d_t)
dm_ar_bg1=lm(d_ar_bg1~1)
acf(dm_ar_bg1$residuals)
dm_ar_bg1$coefficients/sqrt(NeweyWest(dm_ar_bg1,lag=6))

#3-step forecast test
dm_ar_bg3=lm(d_ar_bg3~1)
acf(dm_ar_bg3$residuals)
dm_ar_bg3$coefficients/sqrt(NeweyWest(dm_ar_bg3,lag=6))

#6-step forecast test
dm_ar_bg6=lm(d_ar_bg6~1)
acf(dm_ar_bg6$residuals)
dm_ar_bg6$coefficients/sqrt(NeweyWest(dm_ar_bg6,lag=6))

#12-step forecast test
dm_ar_bg12=lm(d_ar_bg12~1)
acf(dm_ar_bg12$residuals)
dm_ar_bg12$coefficients/sqrt(NeweyWest(dm_ar_bg12,lag=6))


