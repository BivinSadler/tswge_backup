\name{gegenb.wge}
\alias{gegenb.wge}
\title{Calculates Gegenbauer polynomials}
\description{Calculates Gegenbauer polynomials of order n with parameters u and lambda - see (11.9) in Applied Time Series Analysis with R, second edition by Woodward, Gray, and Elliott}
\usage{gegenb.wge(u, d, n)}
\arguments{
  \item{u}{Parameter u in (11.9) Applied Time Series Analysis with R, second edition by Woodward, Gray, and Elliott}
  \item{d}{Parameter lambda in (11.9) Applied Time Series Analysis with R, second edition by Woodward, Gray, and Elliott}
  \item{n}{Order of Gegenbauer polynomial in (11.9)}
}
\details{This function is called by gen.garma.wge}
\value{The coefficients of the nth order Gegenbauer polynomial}
\references{Applied Time Series Analysis with R, second edition by Woodward, Gray, and Elliott}
\author{Wayne Woodward}
\examples{gegenb.wge(u=.8,d=.3,n=6)}
\keyword{ Gegenbauer }

