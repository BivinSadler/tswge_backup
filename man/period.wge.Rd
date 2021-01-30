\name{period.wge}
\alias{period.wge}
\title{Calculate the periodogram}
\description{Given a realization contained in a vector, this function calculates and optionally plots the periodogram in either log or non-log scale}
\usage{period.wge(x, dbcalc = "TRUE", plot = "TRUE")}
\arguments{
  \item{x}{The vector containing the time series realization}
  \item{dbcalc}{if dbcalc=TRUE (default) then the periodogram is calculated in log scale (in dB).  If dbcalc is FALSE then the non-log periodogram is calculated} 
  \item{plot}{if plot=TRUE (default) the periodogram is plotted. If plot=FALSE no plot is created}
}
\value{
  \item{freq}{Frequencies at which the periodogram is calculated}
  \item{pgram}{Periodogram values evaluated at the frequencies in freq}
}
\references{"Applied Time series Analysis with R, 2nd edition" by Woodward, Gray, and Elliott}
\author{Wayne Woodward}
\examples{period.wge(rnorm(100))}

\keyword{ periodogram }
\keyword{ dB }
