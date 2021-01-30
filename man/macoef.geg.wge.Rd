\name{macoef.geg.wge}
\alias{macoef.geg.wge}
\title{Calculate coefficients of the general linear process form of a Gegenbauer process} 
\description{Calculate coefficients of the general linear process form of a Gegenbauer process based on formula (8), page 6 of Ferrara and Guegan(2001).} 
\usage{macoef.geg.wge(u, lambda, trun = 300000)}
\arguments{
  \item{u}{The value of u in the Gegenbauer model}
  \item{lambda}{The value of lambda in the Gegenbauer model}
  \item{trun}{The truncation point of the infinite GLP form}
}
\details{This function is called by gen.geg.wge}
\value{A vector of length trun containing the GLP coefficients}
\references{Applied Time Series Analysis with R, second edition by Woodward, Gray, and Elliott and Ferrara and Guegan(2001)}
\author{Wayne Woodward}
\examples{mageg=macoef.geg.wge(u=.8,lambda=.3)}
\keyword{ Gegenbauer }
\keyword{ GLP }
