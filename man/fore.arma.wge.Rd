\name{fore.arma.wge}
\alias{fore.arma.wge}
\title{Forecast from known model}
\description{Forecasts and associated plots for an ARMA model}
\usage{
fore.arma.wge(x,phi=0,theta=0,n.ahead=5,lastn=FALSE,plot=TRUE,alpha=.05,limits=TRUE)
}
\arguments{
  \item{x}{Realization}
  \item{phi}{AR vector}
  \item{theta}{MA vector}
  \item{n.ahead}{Number of steps ahead}
  \item{lastn}{Logical variable, TRUE means plot forecast for last n.ahead values of realization}
  \item{plot}{Logical variable , TRUE means plot forecasts}
  \item{alpha}{Significance level for prediction limits}
  \item{limits}{Logical variable, TRUE means plot limits}
}
\value{
  \item{f}{Vector of forecasts}
  \item{ll}{Lower limits}
  \item{ul}{Upper limits}
  \item{resid}{Residuals}
  \item{wnv}{White noise variance estimate}
  \item{xbar}{Sample mean of data in x}
  \item{se}{Se for each forecast}
  \item{psi}{psi weights}
  \item{rmse}{RMSE is output if lastn=TRUE}
  \item{mad}{MAD is output if lastn=TRUE}
}
\references{"Applied Time Series Analysis with R, 2nd edition" by Woodward, Gray, and Elliott}
\author{Wayne Woodward}
\examples{data(fig6.1nf)
fore.arma.wge(fig6.1nf,phi=.8,n.ahead=20)}
\keyword{ ARMA }
\keyword{ Forecasts }
