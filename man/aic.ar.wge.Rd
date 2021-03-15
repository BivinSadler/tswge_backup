\name{aic.ar.wge}
\alias{aic.ar.wge}
\title{AR Model Identification for AR models}
\description{AR model identification using either AIC, AICC, or BIC and MLE, Burg or YW}
\usage{aic.ar.wge(x, p = 1:5, type = "aic",method='mle')}
\arguments{
  \item{x}{Realization to be analyzed}
  \item{p}{Range of p values to be considered}
  \item{type}{Type of model identification criterion: aic, aicc, or bic}
  \item{method}{Method used for estimation: MLE, Burg, or YW}
}
\value{
  \item{type}{Criterion used: aic (default), aicc, or bic}
  \item{method}{Estimation method used: MLE, Burg, or YW}
  \item{min_value}{Value of the minimized criterion}
  \item{p}{AR order for selected model}
  \item{phi}{AR parameter estimates for selected model}
  \item{vara}{White noise variance estimate for selected model}
}
\references{"Applied Time Series Analysis with R, 2nd edition" by Woodward, Gray, and Elliott}
\author{Wayne Woodward}
\examples{data(fig3.18a)
          aic.ar.wge(fig3.18a,p=1:5,type='aicc',method='burg')}

\keyword{AIC}
\keyword{AICC}
\keyword{BIC}
\keyword{ML}
\keyword{Burg}
\keyword{YW}
\keyword{Autoregressive}
\keyword{AR Model Identification}
