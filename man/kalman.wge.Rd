\name{kalman.wge}
\alias{kalman.wge}
\title{Kalman filter for simple signal plus noise model}
\description{Kalman filter program to predict, filter, and smooth related to the material in Section 10.6 4 in Applied Time Series Analysis with R, second edition by Woodward, Gray, and Elliott}
\usage{kalman.wge(y, start, gam0, F, gamV, G, gamW)}
\arguments{
  \item{y}{the univariate data set to be analyzed}
  \item{start}{the scalar version of Xo  in item (c) following the state equation (10.47) of the text}
  \item{gam0}{the scalar version of  Gamma(0) discussed in item (c) following the state equation}
  \item{F}{scalar version of the matrix F in the state equation}
  \item{gamV}{the value  Gamma(v) specified in item (b) following the state equation}
  \item{G}{the scalar observation matrix specified in the observation equation as  G(t)}
  \item{gamW}{the variance of the (univariate) white noise denoted by Gamma(w) in item (c) following (10.48)}
}
\value{
  \item{pfs}{a table giving results such as those in Table 10.1 in Woodward, Gray, and Elliott book}
}
\references{Applied Time Series Analysis with R, second edition by Woodward, Gray, and Elliott}
\author{Wayne Woodward}
\note{Requires CRAN package 'astsa'}
\examples{
data(table10.1.signal)
data(table10.1.noise)
spn=table10.1.signal+table10.1.noise
kalman.wge(y=spn,start=0,gam0=1,F=.9,gamV=1,G=1,gamW=.75)
}
\keyword{Kalman filter}
\keyword{predicting}
\keyword{filtering}
\keyword{smoothing}
