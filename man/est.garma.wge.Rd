\name{est.garma.wge}
\alias{est.garma.wge}
\title{Estimate the parameters of a GARMA model.}
\description{This function uses the grid search algorithm discussed in Section 11.5 of Applied Time Series Analysis with R, second edition by Woodward, Gray, and Elliott}
\usage{est.garma.wge(x,low.u,low.lambda,high.u,high.lambda,inc.u,inc.lambda,p.max,nback=500)}
\arguments{
  \item{x}{Realization to be analyzed}
  \item{low.u}{The lower limit for u in the grid search}
  \item{low.lambda}{The lower limit for lambda in the grid search}
  \item{high.u}{The upper limit for u in the grid search}
  \item{high.lambda}{The upper limit for lambda in the grid search}
  \item{inc.u}{The increment, e.g. .01, .001, etc. in the grid search on possible u values}
  \item{inc.lambda}{The increment, e.g. .01, .001, etc. in the grid search on possible lambda values}
  \item{p.max}{Maximum value of p allowed for the AR component of the model}
  \item{nback}{Number of backcasts to be used (see section 11.5 in Applied Time Series Analysis with R, second edition by Woodward, Gray, and Elliott}
}
\details{We assume q=0 and do not allow moving average terms in the model.}
\value{
  \item{u}{Estimate of u}
  \item{lambda}{Estimate of lambda} 
  \item{phi}{Estimates of the pth order AR component of the model where p is some integer from 0 to p.max}
  \item{vara}{The estimated white noise variance}
  \item{aic}{The aic value associated with the final model}
}
\references{Applied Time Series Analysis with R, second edition by Woodward, Gray, and Elliott.  See also Hosking (1984), Gray, Zhang, and Woodward(1989), and Woodward, Cheng, and Gray(1998)}
\author{Wayne Woodward}
\examples{data(llynx)
est.garma.wge(llynx,low.u=.4,high.u=.9,low.lambda=.2,high.lambda=.4,inc.u=.01,inc.lambda=.1,p.max=1)}
\keyword{ GARMA }
\keyword{ Gegenbauer }
\keyword{estimation}

