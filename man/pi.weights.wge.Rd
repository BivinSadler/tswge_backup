\name{pi.weights.wge}
\alias{pi.weights.wge}
\title{Calculate pi weights for an ARMA model}
\description{Given the coefficients of the AR and MA parts of an ARMA model, this function calculates the pi weights }
\usage{pi.weights.wge(phi = 0, theta = 0, lag.max =5)}
\arguments{
  \item{phi}{Vector of AR coefficients (as in "Applied Time Series Analysis with R, 2nd edition" by Woodward, Gray, and Elliott (uses Box and Jenkins notation))}
  \item{theta}{Vector of MA coefficients (as in ATSA and Box Jenkins texts)}
  \item{lag.max}{The function will calculates psi weights pi(1), pi(2), ..., pi(lag.max).  Note that psi(0)=1.}
}
\value{A vector containing pi(1), ..., pi(lag.max)}
\references{"Applied Time Series Analysis with R, 2nd edition" by Woodward, Gray, and Elliott}
\author{Wayne Woodward}
\examples{pi.weights.wge(phi=c(1.2,-.6), theta=.5, lag.max=5)}
\keyword{ Pi weights }

