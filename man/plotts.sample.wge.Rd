\name{plotts.sample.wge}
\alias{plotts.sample.wge}
\title{Plot Data, Sample Autocorrelations, Periodogram, and Parzen Spectral Estimate}
\description{For a given realization, this function plots the data, and calculates and plots the sample autocorrelations, periodogram, and Parzen window spectral estimator in a 2x2 array of plots.}
\usage{plotts.sample.wge(x, lag.max = 25, trunc = 0, arlimits=FALSE,speclimits=c(0,0),
periodogram=FALSE)}
\arguments{
  \item{x}{A vector containing the realization}
  \item{lag.max}{The maximum lag at which to calculate the sample autocorrelations}
  \item{trunc}{The truncation point M for the Parzen spectral estimator. If M=0 theN M=2sqrt(n). If M>0 then M is the value entered}
  \item{arlimits}{Logical variable. TRUE plots 95 percent limit lines on sample autocorrelation plots}
  \item{periodogram}{Logical variable. TRUE plots periodogram, default=FALSE}
  \item{speclimits}{User supplied limits for Parzen spectral density and periodogram, default=function decides limits}
}
  \value{
  \item{xbar}{The sample mean of the realization}
  \item{autplt }{A vector containing sample autocorrelations from 0, 1, ..., aut.lag}
  \item{freq }{A vector containing the frequencies at which the periodogram and window estimate are calculated}
  \item{db }{Periodogram (in dB) calculated at the frequecies in freq}
  \item{freq }{Parzen spectral estimate (in dB) calculated at the frequecies in freq}
}
\references{"Applied Time Series Analysis with R, 2nd edition" by Woodward, Gray, and Elliott}
\author{Wayne Woodward}
\examples{data(wages)
plotts.sample.wge(wages,trunc=0)}
\keyword{ Plot }
\keyword{ Periodogram }
\keyword{ Parzen }