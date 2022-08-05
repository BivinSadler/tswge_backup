\name{expsmooth.wge}
\alias{expsmooth.wge}
\title{Exponential Smoothing}
\description{Performs exponential smoothing on the data in vector x}
\usage{expsmooth.wge(x,alpha=NULL,n.ahead=0,plot=TRUE)}
\arguments{
  \item{x}{Vector containing realization}
  \item{alpha}{Alpha value}
  \item{n.ahead}{Number of steps ahead to forecast}
  \item{plot}{If plot=TRUE then plots of the data along with forecasts}
}
\value{
   \item{alpha}{alpha value used in the smoothing}
   \item{u}{forecasts}
}
\references{"Time Series for Data Science" by Woodward, Sadler, and Robertson}
\author{Wayne Woodward}
\examples{data(wtcrude2020)
expsmooth.wge(wtcrude2020)}
\keyword{ smoothing }
\keyword{ forecasting }
