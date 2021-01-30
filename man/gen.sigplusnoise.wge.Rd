\name{gen.sigplusnoise.wge}
\alias{gen.sigplusnoise.wge}
\title{Generate data from a signal-plus-noise model}
\description{Generate a realization from the model x(t)=coef[1]*cos(2*pi*freq[1]*t+psi[1])+coef[2]*cos(2*pi*freq[2]*t+psi[2])+a(t)}
\usage{gen.sigplusnoise.wge(n,b0,b1=0,coef,freq,psi,phi=0,vara=1,plot=TRUE,sn=0)}
\arguments{
  \item{n}{length of realization to be generated}
  \item{b0}{y intercept of the linear component}
  \item{b1}{slope of the linear component}
  \item{coef}{a 2-component vector specifying the coefficients (if only one cosine term is desired define coef[2]=0)}
  \item{freq}{ a 2-component vector specifying the frequency components (0 to .5)}
  \item{psi}{a 2-component vector specifying the phase shift (0 to 2pi)}
  \item{phi}{a vector of coefficients of the coefficients of the AR noise}
  \item{vara}{vara is the variance of the noise.  NOTE: a(t) is a vector of N(0,WNV) noise generated within the function (default=1)}
  \item{plot}{if TRUE then plot the data generated (default=TRUE)}
  \item{sn}{determines the seed used in the simulation (default=0 indicating new realization each time). sn=positve integer, then the same realization is generated each time}
}
\value{
  \item{x}{realization generated}
}
\references{Applied Time Series Analysis with R, second edition by Woodward, Gray, and Elliott}
\author{Wayne Woodward}
\examples{
x=gen.sigplusnoise.wge(n=100,coef=c(3,1),freq=c(.1,.4),psi=c(0,0),vara=2)
}
\keyword{ signal-plus-noise}
\keyword{ generate}
