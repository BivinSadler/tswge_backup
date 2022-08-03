\name{ample.spec.wge}
\alias{sample.spec.wge}
\title{Smoothed Periodogram using Parzen Window}
\description{This function calculates and optionally plots the smoothed periodogram using the Parzen window. The truncation point may be chosen by the user}
\usage{
sample.spec.wge(x, dbcalc = "TRUE", plot = "TRUE")
}
\arguments{
  \item{x}{Vector containing the time series realization}
  \item{dbcalc}{If dbcalc=TRUE, the calculation is in the log (dB) scale. If FALSE, then non-log calculations are made}
  \item{plot}{If PLOT=TRUE then the smoothed spectral estimate is plotted. If FALSE then no plot is created}
}
\value{
  \item{freq}{The frequencies at which the smoothed periodogram is calculated}
  \item{pzgram }{The smoothed periodogram using the Parzen window}
}
\references{"Applied Time Series Analysis with R, 2nd edition" by Woodward, Gray, and Elliott}
\author{Wayne Woodward}
\examples{sample.spec.wge(rnorm(100))}
\keyword{ Parzen }
\keyword{ Plot }
