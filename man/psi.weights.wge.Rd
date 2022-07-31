\name{psi.weights.wge}
\alias{psi.weights.wge}
\title{Calculate psi weights for an ARMA model}
\description{Given the coefficients of the AR and MA parts of an ARMA model, this function calculates the psi weights }
\usage{psi.weights.wge(phi = 0, theta = 0, lag.max = 5)}
\arguments{
  \item{phi}{Vector of AR coefficients (as in "Applied Time Series Analysis with R, 2nd edition" by Woodward, Gray, and Elliott (uses Box and Jenkins notation))}
  \item{theta}{Vector of MA coefficients (as in ATSA and Box Jenkins texts)}
  \item{lag.max}{The function will calculates psi weights psi(1), psi(2), ..., psi(lag.max).  Note that psi(0)=1.}
}
\value{A vector containing psi(1), ..., psi(lag.max)}
\references{"Applied Time Series Analysis with R, 2nd edition" by Woodward, Gray, and Elliott}
\author{Wayne Woodward}
\examples{psi.weights.wge(phi=c(1.2,-.6), theta=.5, lag.max=5)}
\keyword{ Psi weights }

