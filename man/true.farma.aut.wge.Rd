\name{true.farma.aut.wge}
\alias{true.farma.aut.wge}
\title{True FARMA autocorrelations}
\description{Calculate the autocovariances and autocorrelations and optionally plot the true autocorrlations of a FARMA model}
\usage{
true.farma.aut.wge(d,phi=0,theta=0,lag.max=50,trunc=1000,vara=1,plot=TRUE)
}
\arguments{
  \item{d}{Fractional difference parameter}
  \item{phi}{vector of AR parameters of ARMA part of FARMA model}
  \item{theta}{vector of MA parameters of ARMA part of FARMA model using signs as given in the Woodward, Gray, and Elliott text}
  \item{lag.max}{Maximum lag at which the autocorrelations and autocovariances will be calculated}
  \item{trunc}{Number of terms used in sum}
  \item{vara}{White noise variance}
 \item{plot}{Logical: TRUE=plot, FALSE=no plot}
}
\details{For fractional model use phi=theta=0}
\value{
\item{acf}{Vector of length max.lag+1 containing true autocorrelations at lags 0, 1, ..., lag.max}
\item{acv}{Vector of length max.lag+1 containing true autocovariances at lags 0, 1, ..., lag.max}
}
\references{"Applied Time Series Analysis with R, second editon" by Woodward, Gray, and Elliott}
\author{Wayne Woodward}
\examples{y=true.farma.aut.wge(d=.4,phi=c(0,-.8))}
\keyword{autocorrelations}
\keyword{ fractional}
\keyword{FARMA}
