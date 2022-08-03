\name{slr.wge}
\alias{slr.wge}
\title{Simple Linear Regression}
\description{Uses Base R routine lm to simplify call for SLR where independent variable is automatocally t=1:n}
\usage{slr.wge(x)}
\arguments{
  \item{x}{The TVF data set}
}
\value{
  \item{res}{Residuals}
  \item{b0hat}{Estimate b0 in model y=b0+b1*t+Z}
  \item{b1hat}{Estimate b1}
  \item{pvalue}{pvalue for test:slope=0}
  \item{tstatistic}{tstatistic associated with test:slope=0}
}
\references{Applied Time Series Analysis with R, second edition by Woodward, Gray, and Elliott}
\author{Wayne Woodward}
\examples{x=gen.arma.wge(n=100,phi=.96,sn=10)
                   y=slr.wge(x)}
\keyword{ Regression }
\keyword{ Linear }