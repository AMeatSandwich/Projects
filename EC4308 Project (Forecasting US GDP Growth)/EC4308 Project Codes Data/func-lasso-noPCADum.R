#Here we follow Medeiros et al. (2019) in defining four functions:

#One for forming forecasts using LASSO or Elastic Net model selected on BIC, which will be called
#on each iteration of the rolling window forecasting exercise.

#The second one for producing the series of h-step LASSO forecasts using rolling window.

#The third one will take the resulting coefficients from LASSO, seek out the nonzero ones and rerun OLS with selected predictors to produce post-LASSO forecast.
#The fourth one will call on the third one for producing the series of h-step post-LASSO forecasts using rolling window.

#Inputs for the function:

#1) Data matrix Y: includes all variables

#2) indice - index for dependent variable: 1 for CPI inflation, 2 for PCE inflation

#3) lag - the forecast horizon

#4) alpha - the alpha parameter determining whether we use LASSO (alpha=1) or ElNet (alpha=0.5)

#5) IC - information criterion used to select lambda. Can set to: "bic", "aic" or "aicc"
### define first function: for forming forecasts using lasso selected on BIC, which will
### be called on each iteration of the rolling window forecasting exercise. 
### inputs for this function: 
### 1) data matrix Y, includes all variables 
### 2) indice, index for dependent variable
### 3) lag, the forecast horizon 
### 4) alpha, =1 for lasso 
### 5) ic, information criterion used to select lambda <-  aic, bic, aicc

runlasso=function(Y,indice,lag,alpha,IC){
  aux=embed(Y,4+lag) #create 3 lags + forecast horizon shift (=lag option) for all predictors and INDPRO
  y=aux[,indice] #ICSA variable aligned/adjusted for missing data due to lags
  X=aux[,-c(1:(ncol(Y)*lag))]   # lags of Y (predictors) (lags of ICSA + lags of predictors) corresponding to forecast horizon   
  
  if(lag==1){
    X.out=tail(aux,1)[1:ncol(X)] #retrieve the last  observations of one-step forecast  
  }else{
    X.out=aux[,-c(1:(ncol(Y)*(lag-1)))] #delete first (h-1) columns of aux,
    X.out=tail(X.out,1)[1:ncol(X)] #last observations: y_T,y_t-1...y_t-h
  }
  
  #Here we use the glmnet wrapper written by the authors that does selection on IC:
  model=ic.glmnet(X,y,crit=IC,alpha = alpha) #fit the LASSO (alpha = 1) /ElNet model selected on IC
  
  pred=predict(model,X.out) #generate the forecast 
  
  return(list("model"=model,"pred"=pred)) #return the estimated model and h-step forecast
}


### define second function: this function repeatedly calls the previous function in the rolling window h-step forecasting
### inputs for the function 
### 1) data matrix Y, includes all variables 
### 2) nprev, number of out-of-sample observations (at the end of the sample)
### 3) indice, index for dependent variable -- 1 for ICSA
### 4) lag, the forecast horizon
### 5) alpha, = 1 for lasso
### 6) ic, "bic"/ "aic"/ "aicc"

lasso.rolling.window=function(Y,nprev,indice,lag,alpha,IC){
  save.coef=matrix(NA,nprev,1+4+ncol(Y[,-1])*4 ) #blank matrix for coefficients at each iteration
  save.pred=matrix(NA,nprev,1) #blank for forecasts
  for(i in nprev:1){ #NB: backwards FOR loop: going from 130 down to 1
    Y.window=Y[(1+nprev-i):(nrow(Y)-i),] #define the estimation window (first one: 1 to 491, then 2 to 492 etc.)
    lasso=runlasso(Y.window,indice,lag,alpha,IC) #call the function to fit the LASSO/ElNET selected on IC and generate h-step forecast
    save.coef[(1+nprev-i),]=lasso$model$coef #save estimated coefficients
    save.pred[(1+nprev-i),]=lasso$pred #save the forecast
    cat("iteration",(1+nprev-i),"\n") #display iteration number
  }
  #Some helpful stuff:
  real=Y[,indice] #get actual values
  plot(real,type="l")
  lines(c(rep(NA,length(real)-nprev),save.pred),col="red") #padded with NA for blanks, plot predictions vs. actual
  
  rmse=sqrt(mean((tail(real,nprev)-save.pred)^2)) #compute RMSE
  mae=mean(abs(tail(real,nprev)-save.pred)) #compute MAE (Mean Absolute Error)
  errors=c("rmse"=rmse,"mae"=mae) #stack errors in a vector
  
  return(list("pred"=save.pred,"coef"=save.coef,"errors"=errors)) #return forecasts, history of estimated coefficients, and RMSE and MAE for the period.
}


## Function for doing naive post-LASSO - locate nozero coefficients in inputted LASSO results and rerun OLS ##

#Inputs for the function:

#1) Data matrix Y: includes all variables


#2) indice - index for dependent variable: 1 for CPI inflation, 2 for PCE inflation

#3) lag - the forecast horizon

#4) coef - LASSO coefficients from previous results


################################################################################
### function definitions for post lasso
################################################################################

runpols=function(Y,indice,lag,coef){
  aux=embed(Y,4+lag) #create 12 lags + forecast horizon shift (=lag option)
  y=aux[,indice] #  Y variable aligned/adjusted for missing data due to lags
  X=aux[,-c(1:(ncol(Y)*lag))]   # lags of Y (predictors) corresponding to forecast horizon   
  if(lag==1){
    X.out=tail(aux,1)[1:ncol(X)]  #retrieve last observations if one-step forecast
  }else{
    X.out=aux[,-c(1:(ncol(Y)*(lag-1)))] #delete first (h-1) columns of aux,
    X.out=tail(X.out,1)[1:ncol(X)] #last observations: y_T,y_T-1...y_T-h
  }
  respo=X[,which(coef[-1]!=0)] #find nonzero coefficients in the supplied LASSO coefficient vector and keep corresponding X's
  if(length(respo)==0){ #if no nozero coefficients (full shrinkage)
    model=lm(y ~ 1) #regress on a constant only
  }else{
    model=lm(y ~ respo) #else, run OLS with selected predictors
  }
  coeff=coef(model) #vector of OLS coefficients
  coeff[is.na(coeff)]=0 #if NA set to 0
  pred=c(1,X.out[which(coef[-1]!=0)])%*%coeff #form prediction
  return(list("model"=model,"pred"=pred)) #save OLS model estimates and forecast
}


### function repeatedly calls the post-lasso function in the rolling window h-step forecast 

pols.rolling.window=function(Y,nprev,indice=1,lag,coef){
  
  save.pred=matrix(NA,nprev,1) #blank for forecasts
  for(i in nprev:1){#NB: backwards FOR loop: going from 180 down to 1
    Y.window=Y[(1+nprev-i):(nrow(Y)-i),] #define the estimation window (first one: 1 to 491, then 2 to 492 etc.)
    m=runpols(Y.window,indice,lag,coef[(1+nprev-i),]) #call the function to fit the post-LASSO model and generate h-step forecast
    save.pred[(1+nprev-i),]=m$pred #save the forecast
    cat("iteration",(1+nprev-i),"\n") #display iteration number
  }
  #Some helpful stuff:
  real=Y[,indice] #get actual values
  plot(real,type="l")
  lines(c(rep(NA,length(real)-nprev),save.pred),col="red") #padded with NA for blanks, plot predictions vs. actual
  
  rmse=sqrt(mean((tail(real,nprev)-save.pred)^2)) #compute RMSE
  mae=mean(abs(tail(real,nprev)-save.pred)) #compute MAE (Mean Absolute Error)
  errors=c("rmse"=rmse,"mae"=mae) #stack errors in a vector
  
  return(list("pred"=save.pred,"errors"=errors)) #return forecasts, history of estimated coefficients, and RMSE and MAE for the period.
}

