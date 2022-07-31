\name{ma.pred.wge}
\alias{ma.pred.wge}
\title{Predictive or rolling moving average}
\description{Given a time series in the vector x and order (either an odd or even integer) ma.pred.wge computes a predictive moving average giving 1-step ahead predictions through x(n+1). Optionally, you can specify k-step ahead forecasts beyond the end of the data.}
\usage{ma.pred.wge(x,order=3,n.ahead=1,plot=TRUE)}
\arguments{
  \item{x}{Vector containing original realization}
  \item{order}{Order (odd or even integer) of moving average predictor (default=3)}
  \item{n.ahead}{Number of steps ahead to forecast beyond the end of the data (default=1)}
  \item{plot}{If plot=TRUE then plots of the data and moving average predictors are plotted}
}
\value{
  \item{x}{Original data}
  \item{pred}{Data file showing 1-step ahead predictors up to x(k.ahead)}
  \item{order}{Order (odd or even integer) of the moving average predictor}
  }
\references{"Practical Time Series Analysis with R" by Woodward, Sadler, and Robertson"}
\author{Wayne Woodward}
\examples{data(wtcrude)
sm=ma.pred.wge(x=wtcrude,order=5,n.ahead=10)}
\keyword{ Moving Average Predictor}
\keyword{ Rolling Moving Average}
\keyword{ One-sided Moving Average}
