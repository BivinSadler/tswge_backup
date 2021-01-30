\name{artrans.wge}
\alias{artrans.wge}
\title{Perform Ar transformations}
\description{Given a time series in the vector x, and AR coefs phi1 and phi2, for example, artrans.wge computes y(t)=x(t)-phi1X(t-1)-phi2x(t-2), for t=3, ..., n}
\usage{artrans.wge(x,phi.tr, lag.max=25, plottr = "TRUE")}
\arguments{
  \item{x}{Vector containing original realization}
  \item{phi.tr}{Coefficients of the transformation}
  \item{lag.max}{Max lag (k) for sample autocorrelations}
  \item{plottr}{If plottr=TRUE then plots of the data, transformed data, and sample autocorelations of original and transformed data}
}
\value{Transformed data}
\references{"Applied Time Series Analysis with R, 2nd edition" by Woodward, Gray, and Elliott"}
\author{Wayne Woodward}
\note{For a difference, use phi.tr=1}
\examples{data(wtcrude)
difdata=artrans.wge(wtcrude,phi.tr=1,lag.max=30,plottr=TRUE)}
\keyword{ transformation }
\keyword{ difference }
