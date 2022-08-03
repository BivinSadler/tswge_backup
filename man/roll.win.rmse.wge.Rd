\name{roll.win.rmse.wge}
\alias{roll.win.rmse.wge}
\title{Function to Calculate the Rolling Window RMSE}
\description{This function creates as many "windows" as is possible with the data and calculates an RMSE for each window.  The resulting "rolling window RMSE" is the average of the individual RMSEs from each window.}
\usage{
roll.win.rmse.wge(series, horizon = 2, s = 0, d = 0, phi = 0, theta = 0)
}
\arguments{
  \item{series}{The data}
  \item{horizon}{The number of observations ahead to be forecasted.}
  \item{s}{Order of the seasonal difference, default=1}
  \item{d}{Order of the difference}
  \item{phi}{Vector of AR coefficients}
  \item{theta}{Vector of MA coefficients}
}
\value{
  \item{rwRMSE}{The average of the individual RMSEs of each window}
  \item{numwindows}{The number of windows}
  \item{horizon}{The number of observations ahead to be forecasted.}
  \item{s}{Order of the seasonal difference, default=1}
  \item{d}{Order of the difference}
  \item{phis}{Vector of AR coefficients}
  \item{thetas}{Vector of MA coefficients}
  \item{RMSEs}{Vector of RMSEs ... one for each windwow}
}
\references{"The Time Series Tool Kit"}
\author{Bivin Sadler}

\keyword{ rolling }
\keyword{ window }
\keyword{ RMSE }