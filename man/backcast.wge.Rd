\name{backcast.wge}
\alias{backcast.wge}
\title{Calculate backcast residuals}
\description{This function takes either a fitted (or true) model for the realization x and calculates the residuals using the backcasting procedure}
\usage{backcast.wge(x, phi = 0, theta = 0, n.back = 50)}
\arguments{
  \item{x}{realization}
  \item{phi}{AR coefficients}
  \item{theta}{MA coefficients}
  \item{n.back}{Backcast to X(-n.back)}
}
\value{The n backcast residuals are returned}
\references{Chapter 7 of Applied Time Series Analysis with R, second edition by Woodward, Gray, and Elliott}
\author{Wayne Woodward}
\examples{data(fig6.2nf)
          backcast.wge(fig6.2nf,phi=c(1.2,-.6),theta=.5,n.back=50)}
\keyword{ backcasting }
\keyword{ residuals }
