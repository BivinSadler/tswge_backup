\name{true.garma.aut.wge}
\alias{true.garma.aut.wge}
\title{True GARMA autocorrelations}
\description{Calculate the autocovariances and autocorrelations and optionally plot the true autocorrelations of a 1-factor based on formula(11.25) of "Applied Time Series Analysis with R, second editon" Woodward, Gray, and Elliott}
\usage{true.garma.aut.wge(u,lambda,phi=0,theta=0,lag.max=50,vara=1,plot=TRUE)}
\arguments{
  \item{u}{Parameter u in the GARMA model given in (11.16) of Woodward, Gray, and Elliott text}
  \item{lambda}{Parameter lambda in the GARMA model given in (11.16) of Woodward, Gray, and Elliott text}
  \item{phi}{vector of AR parameters of ARMA part of GARMA model}
  \item{theta}{vector of MA parameters of ARMA part of GARMA model using signs as given in the Woodward, Gray, and Elliott text}
  \item{lag.max}{Maximum lag at which the autocorrelations and autocovariances will be calculated}
  \item{vara}{White noise variance}
 \item{plot}{Logical: TRUE=plot, FALSE=no plot}
}
\details{For Gegenbauer model use phi=theta=0}
\value{
\item{acf}{Vector of length max.lag+1 containing true autocorrelations at lags 0, 1, ..., lag.max}
\item{acv}{Vector of length max.lag+1 containing true autocovariances at lags 0, 1, ..., lag.max}
}
\references{"Applied Time Series Analysis with R, second editon" by Woodward, Gray, and Elliott}
\author{Wayne Woodward}
\examples{y=true.garma.aut.wge(u=.8,lambda=.4,phi=.8)}
\keyword{autocorrelations}
\keyword{ Gegenbauer}
