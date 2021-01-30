\name{parzen.wge}
\alias{parzen.wge}
\title{Smoothed Periodogram using Parzen Window}
\description{This function calculates and optionally plots the smoothed periodogram using the Parzen window. The truncation point may be chosen by the user}
\usage{
parzen.wge(x, dbcalc = "TRUE", trunc = 0, plot = "TRUE")
}
\arguments{
  \item{x}{Vector containing the time series realization}
  \item{dbcalc}{If dbcalc=TRUE, the calculation is in the log (dB) scale. If FALSE, then non-log calculations are made}
  \item{trunc}{if M=0 (default) then the function uses the truncation point 2*sqrt(n). If M>0, then the function uses the given value of M as the truncation point}
  \item{plot}{If PLOT=TRUE then the smoothed spectral estimate is plotted. If FALSE then no plot is created}
}
\value{
  \item{freq}{The frequencies at which the smoothed periodogram is calculated}
  \item{pzgram }{The smoothed periodogram using the Parzen window}
}
\references{"Applied Time Series Analysis with R, 2nd edition" by Woodward, Gray, and Elliott}
\author{Wayne Woodward}
\examples{parzen.wge(rnorm(100))}
\keyword{ Parzen }
\keyword{ Plot }
