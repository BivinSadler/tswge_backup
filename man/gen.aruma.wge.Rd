\name{gen.aruma.wge}
\alias{gen.aruma.wge}
\title{Function to generate an ARUMA (or ARMA or ARIMA) realization}
\description{This function calls arima.sim but an a similar manner to gen.ns.arma.wge and gen.ns.arima.wge but allows for generation of realizations from ARUMA models (see Chapter 5 of "Applied Time Series Analysis with R, 2nd edition" by Woodward, Gray, and Elliott}
\usage{
gen.aruma.wge(n,phi=0,theta=0,d=0,s=0,lambda=0,vara=1,plot=TRUE,sn=0)
}
\arguments{
  \item{n}{Length of realization to be generated}
  \item{phi}{Vector of AR coefficients}
  \item{theta}{Vector of MA coefficients}
  \item{d}{Order of the difference}
  \item{s}{Order of seasonal operator}
  \item{lambda}{Vector of nonstaionary coefficients not associated with d or s (see Def. 5.1(b) in Woodward, Gray, and Elliott text)}
  \item{vara}{White noise variance, default=1}
  \item{plot}{Logical: TRUE=plot, FALSE=no plot}
  \item{sn}{determines the seed used in the simulation. sn=0 produces new/random realization each time. sn=positive integer produces same realization each time}
}
\value{This function generates and (optionally plots) an ARMA or ARIMA or ARUMA realization}
\references{"Applied Time Series Analysis with R, 2nd edition" by Woodward, Gray, and Elliott}
\author{Wayne Woodward}
\examples{gen.aruma.wge(n=100,phi=.7,theta=0, d=1, s=4,lambda=c(1.8,-1),vara=1, plot=TRUE)}
\keyword{ Realization }
\keyword{ ARUMA }
\keyword{ Nonstationary }
