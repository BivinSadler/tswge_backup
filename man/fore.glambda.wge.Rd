\name{fore.glambda.wge}
\alias{fore.glambda.wge}
\title{Forecast using a G(lambda) model}
\description{Find forecasts using a specified G(lambda) model}
\usage{fore.glambda.wge(data.orig,lambda=0,offset=60,phi=0,h=0,n.ahead=10,lastn=TRUE,plot=TRUE)}
\arguments{
  \item{data.orig}{Time series data in the original time scale}
  \item{lambda}{The value of lambda under the Box-Cox time transformation with parameter lambda.}
  \item{offset}{Offset (or shift) value in the G(lambda) model.}
  \item{phi}{Coefficients of the AR component of the AR model fit to the dual data}
  \item{h}{Value of h which will be calculated to produce the desired number of forecasts in the original time scale}
  \item{n.ahead}{Number of values to forecast}
  \item{lastn}{If lastn=TRUE then the last n.ahead values are forecast.  Otherwise, if lastn=FALSE the next n.ahead values are forecast}
  \item{plot}{If plot=TRUE then plots of the data and forecasts are plotted}
}
\details{Forecasts for an AR model fit to the data in the original time scale are also calculated and optionally plotted}
\value{ 
  \item{f.ar}{Forecasts using AR model fit to data in original time}
  \item{f.glam}{Forecasts using AR model fit to the dual and then reinterpolated}
}
\references{Applied Time Series Analysis with R, second edition by Woodward, Gray, and Elliott}
\author{Wayne Woodward}
\examples{data(fig13.2c)
fore.glambda.wge(fig13.2c,lambda=-.4,offset=63,phi=c(0.93,-0.32,-0.15,-0.15,-0.17),n.ahead=30)
}
\keyword{ forecasts }
\keyword{ G(lambda) }
\keyword{ TVF }