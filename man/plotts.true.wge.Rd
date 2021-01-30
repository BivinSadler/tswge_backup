\name{plotts.true.wge}
\alias{plotts.true.wge}
\title{Plot of generated data, true autocorrelations and true spectral density for ARMA model}
\description{For a given ARMA model, this function plots a realization, the true autocorrelations, and the true spectral density. This plot is typical of many plots in Applied Time Series Analysis by Woodward, Gray, and Elliott. For example, see Figure 1.21 and Figure 3.23.}
\usage{
plotts.true.wge(n=100, phi=0, theta=0, lag.max=25, vara = 1)
}
\arguments{
  \item{n}{Length of time series realization to be generated. Default is 100}
  \item{phi}{Vector containing AR parameters}
  \item{theta}{Vector containing MA parameters}
  \item{lag.max}{Maximum lag for calculating and plotting autocorrelations}
  \item{vara}{White noise variance: default=1}
}
\value{
  \item{data }{Realization of length n that is generated from the ARMA model}  
  \item{aut1 }{True autocorrelations from the ARMA model for lags 0 to lag.max}
  \item{acv }{True autocovariances from the ARMA model for lags 0 to lag.max}
  \item{spec}{Spectral density (in dB) for the ARMA model calculated at frequencies f=0, .002, .004, ...., .5}
}
\references{"Applied Time Series Analysis with R, 2nd edition" by Woodward, Gray, and Elliott}
\author{Wayne Woodward}
\note{gvar=g[1], i.e. autocovariance at lag 0}
\examples{plotts.true.wge(n=100, phi=c(1.6,-.9), theta=.8, lag.max=25, vara = 1) }
\keyword{ ARMA }
\keyword{ Autocorrelations }
\keyword{ Spectral density}
