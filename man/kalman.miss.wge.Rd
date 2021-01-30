\name{kalman.miss.wge}
\alias{kalman.miss.wge}
\title{Kalman filter for simple signal plus noise model with missing data}
\description{Kalman function to predict, filter, and smooth in the presence of missing data; see Section 10.6 4 in Applied Time Series Analysis with R} 
\usage{kalman.miss.wge(y,start, gam0, F, gamV, Gtmiss, gamW)}
\arguments{
  \item{y}{the univariate data set to be analyzed}
  \item{start}{the scalar version of X(0) in item (c) following the state equation (10.47) of the text}
  \item{gam0}{the scalar version of  Gamma(0) discussed in item (c) following the state equation}
  \item{F}{scalar version of the matrix F in the state equation}
  \item{gamV}{the value  Gamma(v) specified in item (b) following the state equation}
  \item{Gtmiss}{specifies which items that are missing}
  \item{gamW}{the variance of the (univariate) white noise denoted by Gamma(w) in item (c) following (10.48)}
}
\value{
  \item{pfs}{a table giving results such as those in Table 10.1 in Woodward, Gray, and Elliott book}
}
\references{Applied Time Series Analysis with R, second edition by Woodward, Gray, and Elliott}
\author{Wayne Woodward}
\note{Calls Ksmooth1 in CRAN package 'astsa'}
\examples{
data(table10.1.signal)
data(table10.1.noise)
spn=table10.1.signal+table10.1.noise
n=75
Gtmiss=array(1,dim=c(1,1,n))
Gtmiss[1,1,2]=0
Gtmiss[1,1,5]=0
kalman.miss.wge(y=spn,start=0,gam0=1,F=.9,gamV=1,Gtmiss,gamW=.75)
}
\keyword{Kalman filter}
\keyword{predicting }
\keyword{filtering}
\keyword{smoothing}
\keyword{missing}