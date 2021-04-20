\name{co.wge}
\alias{co.wge}
\title{Cochrane-Orcutt test for trend}
\description{Performs the Cochrane-Orcutt to test for a linear trend in a time series realization.)}
\usage{co.wge(x,maxp=5)}
\arguments{
  \item{x}{Realization}
  \item{maxp}{Maximum AR order allowed for AR model fit to residuals from least squares line}
}
\value{
 \item{z}{ Residuals from the fitted line}
 \item{b0hat}{Estimated y-intercept of the fitted line using the CO method}
 \item{b1hat}{ Estimated slope of the fitted line using the CO method}
 \item{z.order}{Order, p, fit to the residuals}
 \item{z.phi}{Coefficients of the AR model fit to the residuals}
 \item{pvalue}{P-value of the CO test for the significance of the slope}
 \item{tco}{Cochrane-Orcutt test statistic.}  
}
\references{"Applied Time Series Analysis with R, 2nd edition" by Woodward, Gray, and Elliott}
\author{Wayne Woodward}
\examples{data(global.temp)
          co.wge(global.temp,maxp=5)}
\keyword{ Cochrane-Orcutt}
\keyword{ Trend Test }
\keyword{ Regression}
