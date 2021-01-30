\name{gen.garma.wge}
\alias{gen.garma.wge}
\title{Function to generate a GARMA realization}
\description{This function calls gen.geg.wge and arima.sim}
\usage{
gen.garma.wge(n,u,lambda,phi = 0,theta=0,trun=300000,burn_in=600,vara=1,plot=TRUE,sn=0)
}
\arguments{
  \item{n}{the realization length to be generated}
  \item{u}{Parameter u in the GARMA model given in (11.16) of Woodward, Gray, and Elliott text}
  \item{lambda}{Parameter lambda in the GARMA model given in (11.16) of Woodward, Gray, and Elliott text}
  \item{phi}{vector of AR parameters of ARMA part of GARMA model}
  \item{theta}{vector of MA parameters of ARMA part of GARMA model using signs as given ijn the Woodward, Grayu, and Elliott text}
  \item{trun}{the truncation point of the infinite GLP form}
  \item{burn_in}{ is the burning-in period for the simulation}
  \item{vara}{White noise variance, default=1}
  \item{plot}{Logical: TRUE=plot, FALSE=no plot}
  \item{sn}{determines the seed used in the simulation. sn=0 produces new/random realization each time. sn=positive integer produces same realization each time}
}
\value{This function generates and (optionally plots) an GARMA realization}
\references{"Applied Time Series Analysis with R, 2nd edition" by Woodward, Gray, and Elliott}
\author{Wayne Woodward}
\examples{gen.garma.wge(n=100, u=.8,lambda=.4,phi=.9)}
\keyword{ Realization }
\keyword{ GARMA }
\keyword{ Gegenbauer }

