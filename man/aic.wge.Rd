\name{aic.wge}
\alias{aic.wge}
\title{ARMA Model Identification}
\description{ARMA model identification using either AIC, AICC, or BIC}
\usage{aic.wge(x, p = 0:5, q = 0:2, type = "aic")}
\arguments{
  \item{x}{Realization to be analyzed}
  \item{p}{Range of p values to be considered}
  \item{q}{Range of q values to be considered}
  \item{type}{Type of model identification criterion: aic, aicc, or bic}
}
\value{
  \item{type }{Criterion used: aic (default), aicc, or bic}
  \item{min_value }{Value of the minimized criterion}
  \item{p }{AR order for selected model}
  \item{phi }{AR parameter estimates for selected model}
  \item{q }{MA order for selected model}
  \item{theta }{MA parameter estimates for selected model}
  \item{vara }{White noise variance estimate for selected model}
}
\references{"Applied Time Series Analysis with R, 2nd edition" by Woodward, Gray, and Elliott}
\author{Wayne Woodward}
\examples{data(fig3.18a)
          aic.wge(fig3.18a,p=0:5,q=0:1,type='aicc')}

\keyword{ AIC }
\keyword{ AICC }
\keyword{ BIC }
\keyword{ Model Identification}