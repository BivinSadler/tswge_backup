\name{is.glambda.wge}
\alias{is.glambda.wge}
\title{Instantaneous spectrum}
\description{Calculates instantaneous spectrum (in dB) based on a G(lambda) time transformation}
\usage{is.glambda.wge(n, phi = 0, sigma2 = 1, lambda, offset)}
\arguments{
  \item{n}{Length of realization.}
  \item{phi}{Coefficients of AR model fit to dual data.}
  \item{sigma2}{White noise variance}
  \item{lambda}{Lambda in the G(lambda) time transformnation used}
  \item{offset}{Offset in the G(lambda) time transformnation used}
}
\value{Simply a plot of the realization}
\references{"Applied Time Series Analysis with R, 2nd edition" by Woodward, Gray, and Elliott}
\author{Wayne Woodward}
\examples{
is.glambda.wge(n=200,phi=c(.93,-.32,-.15,-.15,-.17),lambda=-.4,offset=63)
}
\keyword{ TVF data }
\keyword{ instantaneous spectrum}
\keyword{time varying spectrum}
\keyword{G(lambda) model}

