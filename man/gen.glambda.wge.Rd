\name{gen.glambda.wge}
\alias{gen.glambda.wge}
\title{Function to generate a g(lambda) realization}
\description{This function generates a g(lambda) TVF realization as discussed in Chapter 13 of Applied Time Series Analysis with R, second edition by Woodward, Gray, and Elliott}
\usage{gen.glambda.wge(n, lambda, phi =0, offset = 20, vara = 1, plot = TRUE, sn = 0)
}
\arguments{
  \item{n}{Length of realization to be generated}
  \item{lambda}{The lambda involved in the g(lambda) time transformation - see Applied Time Series Analysis with R, second edition by Woodward, Gray, and Elliott}  
 \item{phi}{Vector of AR coefficients}
  \item{vara}{White noise variance, default=1}
  \item{offset}{The offset parameter in a g(lambda) process.  See section 13.2 in Applied Time Series Analysis with R, second edition by Woodward, Gray, and Elliott}
  \item{plot}{Logical: TRUE=plot, FALSE=no plot}
  \item{sn}{determines the seed used in the simulation. sn=0 produces new/random realization each time. sn=positive integer produces same realization each time}
}
\value{This function simply generates and (optionally plots) an ARMA realization}
\references{"Applied Time Series Analysis with R, 2nd edition" by Woodward, Gray, and Elliott}
\author{Wayne Woodward}
\examples{gen.glambda.wge(n=500, lambda=0.5,phi=c(1.9,-.99), vara=1, plot=TRUE,sn=0)}
\keyword{ Realization }
\keyword{ g(lambda) process }
\keyword{ TVF }
