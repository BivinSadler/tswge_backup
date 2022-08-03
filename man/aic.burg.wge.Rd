\name{aic.burg.wge}
\alias{aic.burg.wge}
\title{AR Model Identification using Burg Estimates}
\description{AR model identification using either AIC, AICC, or BIC}
\usage{aic.burg.wge(x, p = 1:5, type = "aic")}
\arguments{
  \item{x}{Realization to be analyzed}
  \item{p}{Range of p values to be considered}
  \item{type}{Type of model identification criterion: aic, aicc, or bic}
}
\value{
  \item{type }{Criterion used: aic (default), aicc, or bic}
  \item{min_value }{Value of the minimized criterion}
  \item{p }{AR order for selected model}
  \item{phi }{AR parameter estimates for selected model}
  \item{vara }{White noise variance estimate for selected model}
}
\references{"Applied Time Series Analysis with R, 2nd edition" by Woodward, Gray, and Elliott}
\author{Wayne Woodward}
\examples{data(fig3.18a)
          aic.burg.wge(fig3.18a,p=1:5,type='aicc')}

\keyword{ AIC }
\keyword{ AICC }
\keyword{ BIC }
\keyword{ Model Identification}