\name{ljung.wge}
\alias{ljung.wge}
\title{Ljung-Box Test}
\description{Performs Ljung-Box Test for white noise} 
\usage{ljung.wge(x, K = 24, p = 0, q = 0)}
\arguments{
  \item{x}{Realization to assess for white noise}
  \item{K}{Maximum lag for sample autocorrelations to be used in test}
  \item{p}{If x is a realization of residuals from an ARMA(p,q) fit then p=AR order. Otherwise, p=0}
  \item{q}{If x is a realization of residuals from an ARMA(p,q) fit then q=MA order. Otherwise, q=0}
}
\value{
  \item{test }{Name of test for output: Ljung-Box Test}
  \item{K }{Maximum lag : same as input value}
  \item{chi.square }{Value of chi-square statistic}
  \item{df }{Degrees of freedom = K-p-q}
  \item{pvalue }{pvalue for testing null hypothesis of white noise}
}
\references{"Applied Time Series Analysis with R, 2nd edition" by Woodward, Gray, and Elliott}
\author{Wayne Woodward}
\examples{data(fig1.22a)
ljung.wge(fig1.22a, K=24,p=0,q=0)}
\keyword{ white noise }
\keyword{ Ljung-Box test}
