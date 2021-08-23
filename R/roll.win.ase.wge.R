#Rolling Window ASE Function
# series is the array of the series
# horizon is how far you want to predict into the future
# d is the order of the differencing: (1-B^)^d
# s is the order of the seasonality: (1-B^s)
# phis = order of the stationary AR term
# thetas = order of the invertible MA term

# It simply takes the given horizon and the model in the form of s,d,phis and 
# thetas and figures out how many windows it can create in the data (series) and then calculates the ASE for each window.  
#The output is the average off all the ASEs from each individual window.  

roll.win.ase.wge = function(series, horizon = 1, s = 0, d = 0, phis = 0, thetas = 0)
{
  
  trainingSize = sum(length(phis),length(thetas),s, d) + 1 # sum and plus one is to help backcast.wge, lag.max and ylim plotting issue in fore.arima.wge
  numwindows = length(series)-(trainingSize + horizon) + 1
  RMSEHolder = numeric(numwindows)

  print(paste("Please Hold For a Moment, TSWGE is processing the Rolling Window RMSE with", numwindows, "windows."))
  
  for( i in 1:numwindows)
  {
    
    #invisible(capture.output(forecasts <- fore.arima.wge(series[i:(i+(trainingSize-1))],phi = phis, theta = thetas, s = s, d = d,n.ahead = horizon)))
    forecasts <- fore.arima.wge(series[i:(i+(trainingSize-1))],phi = phis, theta = thetas, s = s, d = d,n.ahead = horizon)
    
    RMSE = sqrt(mean((series[(trainingSize+i):(trainingSize+ i + (horizon) - 1)] - forecasts$f)^2))
    
    RMSEHolder[i] = RMSE
    
  }
  
  RMSEHolder
  hist(RMSEHolder, main = "RMSEs for Individual Windows")
  WindowedRMSE = mean(RMSEHolder)
  
  print("The Summary Statistics for the Rolling Window RMSE Are:")
  print(summary(RMSEHolder))
  print(paste("The Rolling Window RMSE is: ",round(WindowedRMSE,3)))
  invisible(list(rwRMSE = WindowedRMSE, numwindows = numwindows, horizon = horizon, s = s, d = d, phis = phis, thetas = thetas))
}