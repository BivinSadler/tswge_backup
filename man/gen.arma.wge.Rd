\name{gen.arma.wge}
\alias{gen.arma.wge}
\title{Function to generate an ARMA realization}
\description{This function calls arima.sim but with more simple parameter structure for stationary ARMA models}
\usage{
gen.arma.wge(n, phi=0, theta=0, mu=0,vara = 1,plot = TRUE,sn=0)
}
\arguments{
  \item{n}{Length of realization to be generated}
  \item{phi}{Vector of AR coefficients}
  \item{theta}{Vector of MA coefficients}
  \item{vara}{White noise variance, default=1}
  \item{mu}{Theoretical mean, default=0}
  \item{plot}{Logical: TRUE=plot, FALSE=no plot}
  \item{sn}{determines the seed used in the simulation. sn=0 produces new/random realization each time. sn=positive integer produces same realization each time}
}
\value{This function simply generates and (optionally plots) an ARMA realization}
\references{"Applied Time Series Analysis with R, 2nd edition" by Woodward, Gray, and Elliott}
\author{Wayne Woodward}
\examples{gen.arma.wge(n=100, phi=c(1.6,-.9), theta=.8, mu=50,vara=1, plot=TRUE)}

\keyword{ Realization }
\keyword{ ARMA }
