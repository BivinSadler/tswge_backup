\name{pacf.wge}
\alias{pacf.wge}
\title{Compute partial autocorrelations}
\description{Compute partial autocorrelations using either YW (default and the classic al method), Burg, or ML estimates.)}
\usage{pacf.wge(x, lag.max=10 = 2, limits=FALSE, method = "mle")}
\arguments{
  \item{x}{Realization}
  \item{lag.max}{Max AR order}
  \item{limits}{95\% limits}
  \item{method}{Either "mle" (default), "burg", or "yw"}
}
\value{
  \item{method}{Estimation method used: MLE, Burg, or YW}
  \item{pacf}{PACF estimates using estimation method specified}
}
\references{"Time Series for Data Science: Analysis and Forecasting with R" by Woodward,Gray and Elliot}
\author{Wayne Woodward}
\examples{data(sunspot2.0)
          pacf.wge(sunspot2.0,lag.max=10,method='burg')}
\keyword{ PACF }
\keyword{ partial autocorrelations}
\keyword{ Burg }
\keyword{ Yule Walker}
\keyword{ MLE}