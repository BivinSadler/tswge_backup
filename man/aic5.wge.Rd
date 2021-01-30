\name{aic5.wge}
\alias{aic5.wge}
\title{Return top 5 AIC, AICC, or BIC picks}
\description{You may select either AIC, AICC, or BIC to use model identification. Given a range of values for p and q, the program returns the top 5 candidate models.}
\usage{aic5.wge(x, p = 0:5, q = 0:2, type = "aic")}
\arguments{
  \item{x}{Realization to model}
  \item{p}{Range of AR orders to be considered}
  \item{q}{Range of MA orders to be considered}
  \item{type}{Either 'aic' (default, 'aicc', or 'bic')}
}
\value{A list of p,q, and selected criterion for the top 5 models}
\references{"Applied Time Series Analysis with R, 2nd edition" by Woodward, Gray, and Elliott}
\author{Wayne Woodward}
\note{If some model order combinations give explosively nonstationary models, then the program may stop prematurely. You may need to adjust the range of p and q to avoid these models.}
\examples{data(fig3.18a)
          aic5.wge(fig3.18a,p=0:5,q=0:2)}
\keyword{ AIC }
\keyword{ BIC }
\keyword{ AICC}
\keyword{ Model Identification }
