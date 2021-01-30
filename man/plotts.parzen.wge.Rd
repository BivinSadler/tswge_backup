\name{plotts.parzen.wge}
\alias{plotts.parzen.wge}
\title{Calculate and plot the periodogram and Parzen window estimates with differing trunctaion points}
\description{Given a time series contained in the vector x, plotsp.parzen.wge calculates and plots the periodogram and Parzen window estimates at the default truncation point M=2*sqrt(n) and up to 2 additional user specified trunctaion points.}
\usage{plotts.parzen.wge(x, m2=c(0,0))}
\arguments{
  \item{x}{The vector containing the time series realization}
  \item{m2}{A 2-component vector specifying up to 2 additional truncation points}
}
\details{m2=c(10,24) indicates that in addition to the default truncation point, the smoothed spectral estimator
is to be calculated using truncation points 10 and 24, m2=c(0,0) indicates that no additional truncation points are to be used, and m2=c(10,0) indicates the use of one additional truncation point (10)} 
\value{
 \item{freq}{Frequencies at which the periodogram and parzen widow estimates are calculated}
 \item{db}{Periodogram (in dB) calculated at the frequencies in freq}
 \item{dbz}{Parzen window estimate (in dB) calculated at the frequencies in freq using truncation point 2*sqrt(n)}
 \item{dbz1}{Parzen window estimate (in dB) calculated at the frequencies in freq using truncation point m2[1]}
 \item{dbz2}{Parzen window estimate (in dB) calculated at the frequencies in freq using truncation point m2[2]}
}
\references{"Applied Time series Analysis with R, 2nd edition" by Woodward, Gray, and Elliott}
\author{Wayne Woodward}
\examples{data(ss08)
    m2=c(10,50)
    plotts.parzen.wge(ss08,m2)}
\keyword{ periodogram }
\keyword{ window }
\keyword{spectrum}
