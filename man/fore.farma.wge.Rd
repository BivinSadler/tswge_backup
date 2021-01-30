\name{fore.farma.wge}
\alias{fore.farma.wge}
\title{Forecast using a FARMA model}
\description{Find forecasts using a specified FARMA model}
\usage{fore.farma.wge(x, d,  phi, theta = 0, n.ahead = 10, lastn = TRUE, plot = TRUE)}
\arguments{
  \item{x}{Realization  to be analyzed}
  \item{d}{Parameter d in FARMA model}
  \item{phi}{Coefficients of the AR component of the FARMA model}
  \item{theta}{Coefficients of the MA component of the FARMA model}
  \item{n.ahead}{Number of values to forecast}
  \item{lastn}{If lastn=TRUE then the last n.ahead values are forecast.  Otherwise, if lastn=FALSE the next n.ahead values are forecast}
  \item{plot}{If plot=TRUE then plots of the data and forecasts are plotted}
}
\details{Forecasts for an AR model fit to the data are also calculated and optionally plotted}
\value{ 
  \item{ar.fit.order }{Order of the AR model fit to the data}
  \item{ar.fore }{Forecasts based on the AR model}
  \item{farma.fore }{Forecasts based on the FARMA model}
}
\references{Applied Time Series Analysis with R, second edition by Woodward, Gray, and Elliott}
\author{Wayne Woodward}
\examples{fore.farma.wge(Nile, d=.37, phi=0, theta = 0, n.ahead = 30, lastn = TRUE, plot = TRUE)
}
\keyword{ forecasts }
\keyword{ FARMA }
\keyword{ Fractional }
