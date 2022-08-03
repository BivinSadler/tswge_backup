\name{gen.arima.wge}
\alias{gen.arima.wge}
\title{Function to generate an ARIMA (or ARMA) realization}
\description{This function calls arima.sim but with more simple parameter structure for stationary ARIMA (or ARMA) models}
\usage{
gen.arima.wge(n, phi=0, theta=0, d=0,s=0,mu=0,vara=1,plot=TRUE,sn=0)
}
\arguments{
  \item{n}{Length of realization to be generated}
  \item{phi}{Vector of AR coefficients}
  \item{theta}{Vector of MA coefficients}
  \item{d}{Order of the difference}
  \item{s}{Seasonal order}
  \item{vara}{White noise variance, default=1}
  \item{mu}{Theoretical mean of data in x, default=0}
  \item{plot}{Logical: TRUE=plot, FALSE=no plot}
  \item{sn}{determines the seed used in the simulation. sn=0 produces new/random realization each time. sn=positive integer produces same realization each time}
}
\value{This function simply generates and (optionally plots) an ARIMA (or ARMA) realization}
\references{"Applied Time Series Analysis with R, 2nd edition" by Woodward, Gray, and Elliott}
\author{Wayne Woodward}
\examples{gen.arima.wge(n=100, phi=c(1.6,-.9), theta=.8, d=1, vara=1, plot=TRUE)}

\keyword{ Realization }
\keyword{ ARIMA }
\keyword{ Integrated }