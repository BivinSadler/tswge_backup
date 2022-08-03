\name{est.ar.wge}
\alias{est.ar.wge}
\title{Estimate parameters of an AR(p) model}
\description{Estimate parameters of an AR(p) with p assumed known. Outputs residuals (backcast0 and white noise variance estimate.)}
\usage{est.ar.wge(x, p = 2, factor = TRUE, method = "mle")}
\arguments{
  \item{x}{Realization}
  \item{p}{AR order}
  \item{factor}{If TRUE (default) a factor table is printed for the estimated model}
  \item{method}{Either "mle" (default), "burg", or "yw"}
}
\value{
  \item{method}{Estimation method used: MLE, Burg, or YW}
  \item{phi.est}{Estimates of the AR parameters}
  \item{res}{Estimated residuals (using backcasting) based on estimated model}
  \item{avar}{Estimated white noise variance (based on backcast residuals)}
  \item{xbar}{Sample mean of data in x}
  \item{aic}{AIC for estimated model}
  \item{aicc}{AICC for estimated model}
  \item{bic}{BIC for estimated model}
}
\details{The 'type' arument is added for backwards compatabililty and if specified will replace the value specified in the 'method' argument.}
\references{"Applied Time Series Analysis with R, 2nd edition" by Woodward, Gray, and Elliott}
\author{Wayne Woodward}
\examples{data(fig6.1nf)
          est.ar.wge(fig6.1nf,p=1)}
\keyword{ Autoregressive }
\keyword{ Burg }
\keyword{ Yule Walker}
\keyword{ MLE}