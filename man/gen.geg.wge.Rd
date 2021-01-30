\name{gen.geg.wge}
\alias{gen.geg.wge}
\title{Function to generate a Gegenbauer realization}
\description{This function calls macoef.wge}
\usage{
gen.geg.wge(n, u, lambda, trun = 300000,  vara=1 ,sn = 0)
}
\arguments{
  \item{n}{the realization length to be generated}
  \item{u}{Parameter u in the Gegenbauer model given in (11.12) of Woodward, Gray, and Elliott text}
  \item{lambda}{Parameter lambda in the Gegenbauer model given in (11.12) of Woodward, Gray, and Elliott text}
  \item{trun}{the truncation point of the infinite GLP form}
  \item{vara}{White noise variance, default=1}
  \item{sn}{determines the seed used in the simulation. sn=0 produces new/random realization each time. sn=positive integer produces same realization each time}
}
\value{This function generates a Gegenbauer realization}
\details{This function is called by gen.garma.wge and does not have a burn-in time.  Thus, we recommend using est.garma.wge for generating realizations from a Gegenbauer model.}
\references{"Applied Time Series Analysis with R, 2nd edition" by Woodward, Gray, and Elliott}
\author{Wayne Woodward}
\examples{gen.geg.wge(n=100, u=.8,lambda=.4)}
\keyword{ Realization }
\keyword{ Gegenbauer }
\keyword{Long memory }
