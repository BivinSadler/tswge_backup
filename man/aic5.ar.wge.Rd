\name{aic5.ar.wge}
\alias{aic5.ar.wge}
\title{Return top 5 AIC, AICC, or BIC picks for AR model fits}
\description{You may select either AIC, AICC, or BIC to use model identification. You can also used ML, Burg, or Yule-Walker estimates. Given a range of values for p and q, the program returns the top 5 candidate models.}
\usage{aic5.ar.wge(x, p = 0:5, type = "aic",method='mle')}
\arguments{
  \item{x}{Realization to model}
  \item{p}{Range of AR orders to be considered}
  \item{type}{Either 'aic' (default), 'aicc', or 'bic'}
  \item{method}{Either 'MLE' (default), 'Burg', or 'YW'}
}
\value{A list of p, selected criterion for the top 5 models. The identification type and estimation method are printed on the output.}
\references{"Applied Time Series Analysis with R, 2nd edition" by Woodward, Gray, and Elliott}
\author{Wayne Woodward}
\note{If some model order combinations give explosively nonstationary models, then the program may stop prematurely. You may need to adjust the range of p and q to avoid these models.}
\examples{data(fig3.18a)
          aic5.wge(fig3.18a,p=0:5,q=0:2)}
\keyword{ AIC }
\keyword{ BIC }
\keyword{ AICC}
\keyword{ MLE }
\keyword{ Burg }
\keyword{ YW}
\keyword{ Autoregressive}
\keyword{ AR Model Identification }
