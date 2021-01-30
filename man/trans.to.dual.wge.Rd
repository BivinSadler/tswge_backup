\name{trans.to.dual.wge}
\alias{trans.to.dual.wge}
\title{Transforms TVF data set to a dual data set}
\description{Using the specified values for lambda and offset, this function transforms a TVF data set to a dual data set based on a Glambda time transformation.}
\usage{trans.to.dual.wge(x, lambda, offset = 60, h = 0, plot = TRUE)}
\arguments{
  \item{x}{The TVF data set}
  \item{lambda}{The value of lambda in the Glambda time transformation}
  \item{offset}{The value of offset in the Glambda time transformation}
  \item{h}{Scaling variable, initialized at zero, which assures that the dual data set has the same number of points as the original TVF data set} 
  \item{plot}{Logical: TRUE=plot, FALSE=no plot}
}
\value{
  \item{intX}{See intY description}
  \item{intY}{The input realization x is of length n, and the values of x are available at the time points t= 1 to n.  The values intY are n interpolated values of the original time series at the values of intX in the original time scale. The dual data set is obtained by associating the n values of intY with t = 1 to n respectively}
  \item{h}{The output value of the scaling parameter that assures that the dual realization and the original realization are of the same length}
}
\references{Applied Time Series Analysis with R, second edition by Woodward, Gray, and Elliott}
\author{Wayne Woodward}
\examples{data(fig13.2c)
                   y=trans.to.dual.wge(x=fig13.2c,lambda=-.4,offset=63)
}
\keyword{ TVF }
\keyword{ dual }
\keyword{Glambda}
