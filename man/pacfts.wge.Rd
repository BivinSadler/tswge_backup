\name{pacfts.wge}
\alias{pacfts.wge}
\title{Compute partial autocorrelations}
\description{Compute partial autocorrelations using either YW (default and the classical method), Burg, or ML estimates.)}
\usage{pacfts.wge(x,lag.max=5, plot=TRUE,na.action,limits=FALSE,method ='yw')}
\arguments{
  \item{x}{Realization}
  \item{lag.max}{ Max lag}
  \item{plot}{Logical variable}
  \item{na.action}{Not used}
  \item{limits}{Logical variable}
  \item{method}{Either "mle" (default),"burg",or"yw"}
}
\value{
  \item{method}{Estimation method used: MLE, Burg, or YW}
  \item{pacf}{PACF estimates using estimation method specified}
}
\references{"Time Series for Data Science: Analysis and Forecasting with R" by Woodward, Sadler, and Gray}
\author{Wayne Woodward}
\examples{data(sunspot2.0)
          pacfts.wge(sunspot2.0,lag.max=10,method='burg')}
\keyword{ PACF }
\keyword{ partial autocorrelations}
\keyword{ Burg }
\keyword{ Yule Walker}
\keyword{ MLE}