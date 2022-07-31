\name{factor.wge}
\alias{factor.wge}
\title{Produce factor table for a kth order AR or MA model}
\description{This program produces a factor table that reduces a kth order factor into its first and irreducible second order factors as described in Section 3.2.11 of "Applied Time Series Analysis" by Woodward, Gray, and Elliott}
\usage{factor.wge(phi=0, theta=0)}
\arguments{
  \item{phi}{Vector containing the coefficients of the kth order AR factor which is to be factored}
  \item{theta}{Vector containing the coefficients of the kth order MA factor which is to be factored}
}
\value{The only output is the factor table, written by default to the console}
\references{"Applied Time Series Analysis, 2nd edition" by Woodward, Gray, and Elliott}
\author{Wayne Woodward}
\examples{factor.wge(phi=c(-.3,.44,.29,-.378,-.648))}
\keyword{ Factor table }


