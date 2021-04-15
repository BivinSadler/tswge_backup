\name{wbg.boot.wge}
\alias{wbg.boot.wge}
\title{Woodward-Bottone-Gray test for trend}
\description{Performs the Woodward-Bottone-Gray (WBG) bootstrap-based test for a linear trend in a time series realization.)}
\usage{wbg.boot.wge(x,nb=399,alpha=.05,pvalue=TRUE,sn=0)}
\arguments{
  \item{x}{Realization}
  \item{nb}{The number of Bootstrap replications (default is 399)}
  \item{alpha}{The significance level of the test (default is .05)}
  \item{pvalue}{Logical variable. TRUE(default) prints out the p-value of the test.}
  \item{sn}{Sets the seed for the simulations (default = 0))}
}
\value{
  \item{p}{AR order used for the bootstrap simulations}
  \item{phi}{The AR coefficients of the AR model fit to data}
  \item{pv}{The p-value of the test}
}
\references{"Applied Time Series Analysis with R, 2nd edition" by Woodward, Gray, and Elliott}
\author{Wayne Woodward}
\examples{data(global.temp)
          wbg.boot.wge(global.temp)}
\keyword{ Woodward-Bottone-Gray}
\keyword{ Trend Test }
\keyword{ Regression}
\keyword{Bootstrap}


