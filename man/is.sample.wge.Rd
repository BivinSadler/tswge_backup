\name{is.sample.wge}
\alias{is.sample.wge}
\title{Sample instantaneous spectrum based on periodogram}
\description{Calculates sample instantaneous spectrum (in dB) based on a G(lambda) time transformation}
\usage{is.sample.wge(data,  lambda, offset)}
\arguments{
  \item{data}{Realization to be analyzed.}
  \item{lambda}{Lambda in the G(lambda) time transformation used}
  \item{offset}{Offset in the G(lambda) time transformation used}
}
\value{Simply a plot of the realization}
\references{"Applied Time Series Analysis with R, 2nd edition" by Woodward, Gray, and Elliott}
\author{Wayne Woodward}
\examples{
data(ss08)
is.sample.wge(data=ss08,lambda=-.4,offset=63)
}
\keyword{ TVF data }
\keyword{ instantaneous spectrum}
\keyword{time varying spectrum}
\keyword{G(lambda) model}