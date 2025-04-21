#Here we follow Medeiros et al. (2019) in defining two functions:

#One for forming forecasts using AR(p) model selected on BIC, which will be called
#on each iteration of the rolling window forecasting exercise.

#The other one for producing the series of h-step forecasts using rolling window.

#Inputs for the function:

#1) Data matrix Y: includes all variables

#2) indice - index for dependent variable: 1 for CPI inflation, 2 for PCE inflation

#3) lag - the forecast horizon

################################################################################
### function definitions for random forest
################################################################################

### this function runs rf on one rolling window

runrf=function(Y,indice,lag){
  
  aux=embed(Y,4+lag) #create 4 lags + forecast horizon shift (=lag option)
  y=aux[,indice] #  Y variable aligned/adjusted for missing data due do lags
  X=aux[,-c(1:(ncol(Y)*lag))]  # lags of Y (predictors) corresponding to forecast horizon 
  if(lag==1){
    X.out=tail(aux,1)[1:ncol(X)]   #retrieve the last  observations if one-step forecast
  }else{
    X.out=aux[,-c(1:(ncol(Y)*(lag-1)))] #delete first (h-1) columns of aux,
    X.out=tail(X.out,1)[1:ncol(X)]  #last observations: y_T,y_t-1...y_t-h
  }
  model=randomForest(X,y,importance = TRUE) #fit the random forest on default settings
  pred=predict(model,X.out) #generate forecast
  return(list("model"=model,"pred"=pred)) #return the estimated model and h-step forecast
}

### this function repeatedly calls the previous function in the rolling window h-step forecasting
rf.rolling.window=function(Y,nprev,indice,lag){
  save.importance=list() #blank for saving variable importance
  save.pred=matrix(NA,nprev,1) ##blank for forecasts
  for(i in nprev:1){#NB: backwards FOR loop: going from 130 down to 1
    Y.window=Y[(1+nprev-i):(nrow(Y)-i),] #define the estimation window (first one: 1 to 392, then 2 to 393 etc.)
    rf=runrf(Y.window,indice,lag)#call the function to fit the Random Forest and generate h-step forecast
    save.pred[(1+nprev-i),]=rf$pred #save the forecast
    save.importance[[i]]=importance(rf$model) #save variable importance
    cat("iteration",(1+nprev-i),"\n") #display iteration number
  }
  #Some helpful stuff:
  real=Y[,indice]#get actual values
  plot(real,type="l")
  lines(c(rep(NA,length(real)-nprev),save.pred),col="red") #padded with NA for blanks, plot predictions vs. actual
  rmse=sqrt(mean((tail(real,nprev)-save.pred)^2)) #compute RMSE
  mae=mean(abs(tail(real,nprev)-save.pred)) #compute MAE (Mean Absolute Error)
  errors=c("rmse"=rmse,"mae"=mae) #stack errors in a vector
  return(list("pred"=save.pred,"errors"=errors,"save.importance"=save.importance)) #return forecasts, history of variable importance, and RMSE and MAE for the period.
}

