\name{est.ar.wge}
\alias{est.ar.wge}
\title{Estimate parameters of an AR(p) model}
\description{Estimate parameters of an AR(p) with p assumed known. Outputs residuals (backcast0 and white noise variance estimate.)}
\usage{est.ar.wge(x, p = 2, factor = TRUE, type = "mle")}
\arguments{
  \item{x}{Realization}
  \item{p}{AR order}
  \item{factor}{If TRUE (default) a factor table is printed for the estimated model}
  \item{type}{Either "burg" (default), "yw", or "mle"}
}
\value{
  \item{phi.est}{Estimates of the AR parameters}
  \item{res}{Estimated residuals (using backcasting) based on estimated model}
  \item{vara}{Estimated white noise variance (based on backcast residuals)}
  \item{aic}{AIC for estimated model}
  \item{aicc}{AICC for estimated model}
  \item{bic}{BIC for estimated model}
}
\references{"Applied Time Series Analysis with R, 2nd edition" by Woodward, Gray, and Elliott}
\author{Wayne Woodward}
\examples{data(fig6.1nf)
          est.ar.wge(fig6.1nf,p=1)}
\keyword{ Autoregressive }
\keyword{ Burg }
\keyword{ Yule Walker}
\keyword{ MLE}