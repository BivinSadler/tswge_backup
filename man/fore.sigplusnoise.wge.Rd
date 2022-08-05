\name{fore.sigplusnoise.wge}
\alias{fore.sigplusnoise.wge}
\title{Forecasting signal plus noise models}
\description{Forecast models of the form line plus AR noise or cosine plus AR noise with known frequency}
\usage{fore.sigplusnoise.wge(x,linear=TRUE,method="mle",freq=0,max.p=5,
n.ahead=10,lastn=FALSE,plot=TRUE,alpha=.05,limits=TRUE)}
\arguments{
  \item{x}{The variable containing the realization to be analyzed}
  \item{linear}{If TRUE then the program forecasts a line plus noise model.  If FALSE the model is cosine plus noise}
  \item{method}{Estimation method}
  \item{freq}{Frequency of the cosine term. freq is ignored when using line plus noise}
  \item{max.p}{Max value of p for the ARp model fit to the noise}
  \item{n.ahead}{The number of steps ahead to forecast}
  \item{lastn}{If TRUE then the function forecasts the last n.ahead values of the realization. If FALSE the the forecasts are for n.ahead steps beyond the end of the realization}
  \item{plot}{If TRUE then the forecasts and realization are plotted}
  \item{alpha}{Significance level}
  \item{limits}{If TRUE the forecast limits calculated and plotted}
}
\value{
  \item{f}{The n.ahead forecasts}
  \item{ll}{The lower limits for the forecasts. zeros are returned if limits were not requested}
  \item{ul}{The upper limits for the forecasts. zeros are returned if limits were not requested}
  \item{res}{Residuals}
  \item{wnv}{The estimated white noise variance based on the residuals}
  \item{se}{se is the estimated standard error of the k step ahead forecast. zeros are returned if limits were not requested}
  \item{xi}{xi is the kth psi weight associated with the fitted AR model and used to calculate the se above.  Note that psi0 is1. zeros are returned if limits were not requested}
}
\references{"Applied Time Series Analysis with R, 2nd edition" by Woodward, Gray, and Elliott}
\author{Wayne Woodward}
\examples{
data(llynx)
llynx.for=fore.sigplusnoise.wge(llynx,linear=FALSE,freq=.1,max.p=5,n.ahead=20)
}
\keyword{ signal plus noise }
\keyword{ forecast }
