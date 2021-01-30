\name{trans.to.original.wge}
\alias{trans.to.original.wge}
\title{Transforms  dual data set back to original time scale}
\description{Using the specified values for lambda and offset, this function transforms a dual data set, based on a Glambda time transformation, back to the original time scale}
\usage{trans.to.original.wge(xd, lambda, offset, h, plot = TRUE)}
\arguments{
  \item{xd}{The dual data set}
  \item{lambda}{The value of lambda in the Glambda time transformation}
  \item{offset}{The value of offset in the Glambda time transformation}
  \item{h}{Scaling variable obtained as output from transform.to.dual.wge that assures that the dual data set has the same number of points as the origuinal TVF data set} 
  \item{plot}{Logical: TRUE=plot, FALSE=no plot}
}
\value{Returns the y values to be plotted at time points  t=1 to n that approximate the original TVF data set}
\references{Applied Time Series Analysis with R, second edition by Woodward, Gray, and Elliott}
\author{Wayne Woodward}
\examples{data(fig13.2c)
                   yd=trans.to.dual.wge(fig13.2c,lambda=-.4,offset=63)
                   yo=trans.to.original.wge(yd$intY,lambda=-.4,offset=63,h=yd$h)
}
\keyword{ TVF }
\keyword{ dual }
\keyword{Glambda}
