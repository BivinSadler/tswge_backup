\name{est.glambda.wge}
\alias{est.glambda.wge}
\title{Estimate the value of lambda and offset to produce a stationary dual.}
\description{This function uses the technique discussed in Section 13.3.3 of Applied Time Series Analysis with R, second edition by Woodward, Gray, and Elliott to find the g(lambda) time transformation that most nearly transforms the data to a stationary dual.}
\usage{est.glambda.wge(data, lambda.range = c(0, 1), offset.range = c(0, 100))}
\arguments{
  \item{data}{Vector containing the TVF realization to be analyzed}
  \item{lambda.range}{Range of lambda values considered in the search}
  \item{offset.range}{Range of offset values considered in the search}
}
\value{
  \item{Q}{A listing of lambda values within the range and offsets for each lambda that provided the best dual.  Also a listing of the test statistic, Q, to be minimized}
  \item{best.lambda}{See description of best.offset below}
  \item{best.offset}{best.lambda and best.offset are the lambda-offset pair that produced the most stationary dual according to the Q criterion}
}
\references{Applied Time Series Analysis with R, second edition by Woodward, Gray, and Elliott and Jiang, Gray, and Woodward(2006)}
\author{Wayne Woodward}
\examples{data(ss08)
                   est.glambda.wge(ss08,lambda.range=c(-1,1),offset.range=c(0,100))
}
\keyword{ TVF }
\keyword{ g(lambda)}
