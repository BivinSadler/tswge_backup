\name{true.arma.aut.wge}
\alias{true.arma.aut.wge}
\title{True ARMA autocorrelations}
\description{R function to calculate the autocovariances and autocorrelations and optionally plot the true autocorrelations of a stationary ARMA model}
\usage{
true.arma.aut.wge(phi = 0, theta = 0, lag.max = 25, vara = 1,plot=TRUE)
}
\arguments{ 
  \item{phi}{Vector containing AR coefficients}
  \item{theta}{Vector containing MA coefficients}
  \item{lag.max}{Maximum lag at which to calculate the true autocorrelations}
  \item{vara}{White noise variance of the ARMA model}
 \item{plot}{Logical: TRUE=plot, FALSE=no plot}
}
\value{
\item{acf}{Vector of length max.lag+1 containing true autocorrelations at lags 0, 1, ..., lag.max}
\item{acv}{Vector of length max.lag+1 containing true autocovariances at lags 0, 1, ..., lag.max}
}
\references{"Applied Time Series Analysis with R, 2nd edition" by Woodward, Gray, and Elliott}
\author{Wayne Woodward}
\examples{true.arma.aut.wge(phi=c(1.6,-.9),theta=-.8,lag.max=15,vara=1) }
\keyword{ Autocorrelations }
\keyword{ ARMA }
