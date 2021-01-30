\name{true.arma.spec.wge}
\alias{true.arma.spec.wge}
\title{True ARMA Spectral Density}
\description{R function to calculate and optionally plot the spectral density of a stationary ARMA model}
\usage{true.arma.spec.wge(phi=0,theta=0, vara=1,plot=TRUE)}

\arguments{ 
  \item{phi}{Vector containing AR coefficients}
  \item{theta}{Vector containing MA coefficients}
  \item{vara}{White noise variance of the ARMA model}
  \item{plot}{Logical: TRUE=plot, FALSE=no plot}
}
\value{
 \item{f}{Frequencies at which true spectral density is evaluated: 0, 1/500, 2/500, ..., .5}
  \item{spec }{True spectral density calculated at the frequencies in f}
}
\references{"Applied Time Series Analysis with R, 2nd edition" by Woodward, Gray, and Elliott}
\author{Wayne Woodward}
\examples{true.arma.spec.wge(phi=c(1.6,-.9), theta=.7) }
\keyword{ spectral density }
\keyword{ ARMA }
