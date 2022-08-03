\name{fore.arima.wge}
\alias{fore.arima.wge}
\title{Function for forecasting from known model which may have (1-B)^d and/or seasonal factors}
\description{This function calculates forecasts from a known model that may have stationary ARMA components as well as (1-B)^dand/or seasonal factors}
\usage{
fore.arima.wge(x,phi=0,theta=0,d=0,s=0,n.ahead=5,lastn=FALSE,plot=TRUE,alpha=.05,limits)
}
\arguments{
  \item{x}{Realization to be forecast from}
  \item{phi}{Vector containing stationary AR parameters}
  \item{theta}{Vector containing MA parameters}
  \item{d}{Order of difference}
  \item{s}{Seasonal order}
  \item{n.ahead}{Number of steps ahead to forecast}
  \item{lastn}{Logical, lastn=TRUE plots forecasts for the last n.ahead values in the realization}
  \item{plot}{Logical, plot=TRUE plots forecasts}
  \item{alpha}{Significance level for prediction limits}
  \item{limits}{Logical, limits=TRUE plots prediction limits}
}
\value{
  \item{f}{Vector of forecasts}
  \item{ll}{Lower limits}
  \item{ul}{Upper limits}
  \item{resid}{Residuals}
  \item{wnv}{White noise variance estimate}
  \item{xbar}{Sample mean of data in x}
  \item{se}{Se for each forecast}
  \item{psi}{Psi weights}
  \item{ptot}{Total order of all AR components, phi, d, and s}
  \item{phtot}{Coefficients after multiplying all stationary and nonstationary coponents on the AR side of the equation}
}
\references{"Applied Time Series Analysis with R, 2nd edition" by Woodward, Gray, and Elliott}
\author{Wayne Woodward}
\examples{data(airline)
          x=log(airline)
          phi12=c(-.36,-.05,-.14,-.11,.04,.09,-.02,.02,.17,.03,-.1,-.38)
          s=12
          d=1
          fore.arima.wge(x,phi=phi12,d=1,s=12,n.ahead=12,limits=FALSE)}
\keyword{ Forecasts }
\keyword{ ARIMA}
\keyword{ Seasonal}
