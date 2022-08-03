\name{roll.win.rmse.nn.wge}
\alias{roll.win.rmse.nn.wge}
\title{Function to Calculate the Rolling Window RMSE}
\description{This function creates as many "windows" as is possible with the data and calculates an RMSE for each window.  The resulting "rolling window RMSE" is the average of the individual RMSEs from each window.}
\usage{
roll.win.rmse.nn.wge(series, horizon = 1, fit_model)
}
\arguments{
  \item{series}{The data}
  \item{horizon}{The number of observations ahead to be forecasted.}
  \item{fit_model}{The mlp object (model) to be evaluated.  This model will have been fit before the call to this function.}

}
\value{
  \item{rwRMSE}{The average of the individual RMSEs of each window}
  \item{numwindows}{The number of windows}
  \item{horizon}{The number of observations ahead to be forecasted.}
}

\references{"The Time Series Tool Kit"}
\author{Bivin Sadler}

\keyword{ rolling }
\keyword{ window }
\keyword{ RMSE }