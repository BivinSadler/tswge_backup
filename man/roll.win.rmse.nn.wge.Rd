\name{roll.win.rmse.nn.wge}
\alias{roll.win.rmse.nn.wge}
\title{Function to Calculate the Rolling Window RMSE}
\description{This function creates as many "windows" as is possible with the data and calculates an ASE for each window.  The resulting "rolling window RMSE" is the average of the individual RMSEs from each window.}
\usage{
roll.win.rmse.nn.wge = function(series, horizon = 1, s = 0, d = 0, phi = 0, theta = 0)
}
\arguments{
  \item{series}{The data}
  \item{horizon}{The number of observations ahead to be forecasted.}

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