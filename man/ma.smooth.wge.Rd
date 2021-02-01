\name{ma.smooth.wge}
\alias{ma.smooth.wge}
\title{Centered Moving Average Smoother}
\description{Given a time series in the vector x and order (either an odd or even integer) ma.smooth.wge computes a centered moving average smoother and optionally plots the data and smoothed data}
\usage{ma.smooth.wge(x,order=3,plot=TRUE)}
\arguments{
  \item{x}{Vector containing original realization}
  \item{order}{Order (odd or even integer) of moving average smoother}
  \item{plot}{If plot=TRUE then plots of the data and smoothed data are plotted}
}
\value{
  \item{x}{Original data}
  \item{smooth}{Data after application of centered average filter.l}
  \item{order}{Order (odd or even integer) of the smoother}
  }
\references{"Practical Time Series Analysis with R" by Woodward, Sadler, and Robertson"}
\author{Wayne Woodward}
\examples{data(wtcrude)
sm=ma.smooth.wge(x=wtcrude,order=5)}
\keyword{ Moving Average Smoother}
\keyword{ Centered}
