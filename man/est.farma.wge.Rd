\name{est.farma.wge}
\alias{est.farma.wge}
\title{Estimate the parameters of a FARMA model.}
\description{This function uses the grid search algorithm discussed in Section 11.5 of Applied Time Series Analysis with R, second edition by Woodward, Gray, and Elliott}
\usage{est.farma.wge(x, low.d, high.d, inc.d, p.max, nback = 500)}
\arguments{
  \item{x}{Realization to be analyzed}
  \item{low.d}{The lower limit for d in the grid search}
  \item{high.d}{The upper limit for d in the grid search}
  \item{inc.d}{The increment, e.g. .01, .001, etc. in the grid search}
  \item{p.max}{Maximum value of p allowed for the AR component of the model}
  \item{nback}{Number of backcasts to be used (see section 11.5 in Applied Time Series Analysis with R, second edition by Woodward, Gray, and Elliott}
}
\details{We assume q=0 and do not allow moving average terms in the model.}
\value{
  \item{d}{Estimate of d}
  \item{phi}{Estimates of the pth order AR component of the model where p is some integer from 0 to p.max}
  \item{vara}{The estimnated white noise variance}
  \item{aic}{The aic value associated with the final model}
}
\references{Applied Time Series Analysis with R, second edition by Woodward, Gray, and Elliott.  See also Hosking (1984)}
\author{Wayne Woodward}
\examples{est.farma.wge(Nile,low.d=.1,high.d=.5,inc.d=.01,p.max=3)}
\keyword{ FARMA }
\keyword{ estimation }
