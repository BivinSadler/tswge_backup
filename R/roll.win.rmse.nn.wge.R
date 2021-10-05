#  big problem with passing mlp an old model is that it doesn't always change the difforder.
# Sometimes it overrides the old model difforder so there is code to manually force difforder to be 0
# by overwriting the instantiated mlp object (in the loop) wiht the old model difforder.  
# There is an excel file (RollWindowNNCheck.xlsx) that does the MSE calculation by hand to double check but the below 
# code compares the mlp function MSE to the rolling window MSE for horizon 1 which should be the MSE.  
# Two issues arise
# 1. The code automatically checks for trend and uses the forecast::ets function which requires a series of at least 8.
# Therefore, the rolling window RMSE won't match exactly because it is always missing the first 8 windows. (But it should be close).
# 2. The difforder is not always taken from the old model and windows with perceived trend will difference 
# when there shouldn't be one. The fix for this is to manually force difforder to
# overwrite the instantiated mlp object's difforder (in the loop) with the old model difforder.


roll.win.rmse.nn.wge = function(series, horizon = 1, fit_model)
{

if(fit_model$sdummy == FALSE)
{
  lags = fit_model$lags
  if(is.null(fit_model$difforder)){d = 0}
  else{
  d = max(fit_model$difforder)
  }
  s = frequency(series)
  trainingSize = max(lags+sum(d)+2,9) #Long story but forecast::ets does not like series smaller than 8 unit and nnfor calls "trendcheck" when difforder = 0 which returns error for lenth < 8.
  print(trainingSize)
  numwindows = length(series)-(trainingSize + horizon) + 1
  print(numwindows)
  RMSEHolder = numeric(numwindows)
  MSEHolder = numeric(numwindows)
  #ForecastHolder = numeric(numwindows)
  
  print(paste("Please Hold For a Moment, TSWGE is processing the Rolling Window RMSE with", numwindows, "windows."))
  
  for( i in 1:numwindows)
  {
    fit = mlp(ts(series[i:(i+(trainingSize-1))]), model = fit_model)
    fit$minmax = fit_model$minmax
    fit$difforder = fit_model$difforder
    forecasts = predict(fit,h = horizon)$mean
    MSE = mean((series[(trainingSize+i):(trainingSize+ i + (horizon) - 1)] - forecasts)^2)
    MSEHolder[i] = MSE
    #ForecastHolder[i] = forecasts
  }
  
  RMSEHolder = sqrt(MSEHolder)
  hist(RMSEHolder, main = "RMSEs for Individual Windows")
  WindowedMSE = mean(MSEHolder)
  WindowedRMSE = mean(RMSEHolder)
  
  print("The Summary Statistics for the Rolling Window RMSE Are:")
  print(summary(RMSEHolder))
  print(paste("The Rolling Window RMSE is: ",round(WindowedRMSE,3)))
  invisible(list(rwRMSE = WindowedRMSE, rwMSE = WindowedMSE, numwindows = numwindows, horizon = horizon, s = s, d = d))
} # end big if
else   # This is if there is seasonality ... the problem is that if you give the model, mlp() still recalculates the seasonal dummies.
{
    lags = fit_model$lags
    if(is.null(fit_model$difforder)){d = 0}
    else{
      d = max(fit_model$difforder)
    }
    s = frequency(series)
    trainingSize = 20 # It will recalculate seasonal dummies... need large portion of data
    print(trainingSize)
    #numwindows = length(series)-(trainingSize + horizon) + 1
    numwindows = 1
    print(numwindows)
    RMSEHolder = numeric(numwindows)
    MSEHolder = numeric(numwindows)
    #ForecastHolder = numeric(numwindows)
    
    print(paste("Please Hold For a Moment, TSWGE is processing the Rolling Window RMSE with", numwindows, "windows."))
    print(paste("Seasonal MLP models will be evaluated only on the last window."))
    
    for( i in 1:numwindows)
    {
      fit = mlp(ts(series[1:(length(series)-horizon)], frequency = s), model = fit_model)
      fit$minmax = fit_model$minmax
      fit$difforder = fit_model$difforder
      #print(fit)
      #print(fit$ff.det)
      forecasts = predict(fit,h = horizon)$mean
      MSE = mean((series[(length(series)-(horizon)+1):length(series)] - forecasts)^2)
      MSEHolder[i] = MSE
      #ForecastHolder[i] = forecasts
    }
    
    RMSEHolder = sqrt(MSEHolder)
    hist(RMSEHolder, main = "RMSEs for Individual Windows")
    WindowedMSE = mean(MSEHolder)
    WindowedRMSE = mean(RMSEHolder)
    
    print("The Summary Statistics for the Rolling Window RMSE Are:")
    print(summary(RMSEHolder))
    print(paste("The Rolling Window RMSE is: ",round(WindowedRMSE,3)))
    invisible(list(rwRMSE = WindowedRMSE, rwMSE = WindowedMSE, numwindows = numwindows, horizon = horizon, s = s, d = d, MH= MSEHolder))
  } # end big else
} # end function 
