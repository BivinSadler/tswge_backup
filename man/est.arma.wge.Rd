\name{est.arma.wge}
\alias{est.arma.wge}
\title{Function to calculate ML estimates of parameters of stationary ARMA models}
\description{This function calculates ML estimates, computes residuals (using backcasting), estimates white noise variance for a stationary ARMA model}
\usage{est.arma.wge(x, p = 0, q = 0, factor = TRUE)}
\arguments{
  \item{x}{The realization.}
  \item{p}{The autoregressive order}
  \item{q}{the moving average order}
  \item{factor}{Logical variable.  factor=TRUE (default) plots a factor table for estimated AR-part of model}
}
\details{This function uses arima from base SAS and is written similarly to itsmr function arma}
\value{
  \item{phi}{ML estimates of autoregressive parameters}
  \item{theta}{ML estimates of moving average parameters}
  \item{res}{Residuals (calculated using backcasting)}
  \item{avar}{Estimate of white noise variance based on backcast residuals}
  \item{se.phi}{Standard errors of the AR parameter estimates}
  \item{se.theta}{Standard errors of the MA parameter estimates}
  \item{aic}{AIC for estimated model}
  \item{aicc}{AICC for estimated model}
  \item{bic}{BIC for estimated model}
}
\references{"Applied Time Series Analysis with R, 2nd edition" by Woodward, Gray, and Elliott}
\author{Wayne Woodward}
\note{Requires CRAN package 'itsmr'. The program is based on arima from base R and arma from 'itsmr'}
\examples{data(fig6.2nf)
          est.arma.wge(fig6.2nf,p=2,q=1)}
\keyword{ Maximum likelihood }
\keyword{ Backcasting }