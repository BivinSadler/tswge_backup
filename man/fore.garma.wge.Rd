\name{fore.garma.wge}
\alias{fore.garma.wge}
\title{Forecast using a GARMA model}
\description{Find forecasts using a specified GARMA model}
\usage{fore.garma.wge(x,u,lambda,phi,theta=0,n.ahead=10,lastn=TRUE,plot=TRUE)}
\arguments{
  \item{x}{Realization  to be analyzed}
  \item{u}{Parameter u in GARMA model}
  \item{lambda}{Parameter lambda in GARMA model}
  \item{phi}{Coefficients of the AR component of the GARMA model}
  \item{theta}{Coefficients of the MA component of the GARMA model}
  \item{n.ahead}{Number of values to forecast}
  \item{lastn}{If lastn=TRUE then the last n.ahead values are forecast.  Otherwise, if lastn=FALSE the next n.ahead values are forecast}
  \item{plot}{If plot=TRUE then plots of the data and forecasts are plotted}
}
\details{Forecasts for an AR model fit to the data are also calculated and optionally plotted}
\value{ 
  \item{ar.fit.order }{Order of the AR model fit to the data}
  \item{ar.fore }{Forecasts based on the AR model}
  \item{garma.fore }{Forecasts based on the GARMA model}
}
\references{Applied Time Series Analysis with R, second edition by Woodward, Gray, and Elliott}
\author{Wayne Woodward}
\examples{data(llynx)
fore.garma.wge(llynx,u=.796,lambda=.4,phi=.51,theta=0,n.ahead=30,lastn=TRUE,plot=TRUE)
}
\keyword{ forecasts }
\keyword{ GARMA }
\keyword{ Gegenbauer }
