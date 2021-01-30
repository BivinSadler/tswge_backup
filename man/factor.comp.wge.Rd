\name{factor.comp.wge}
\alias{factor.comp.wge}
\title{Create a factor table and AR components for an AR realization}
\description{This program finds the ML estimates of a specified order, then prints a factor table for the estimated model and prints and plots the additive components}
\usage{
factor.comp.wge(x, aic = FALSE, p, ncomp)
}
\arguments{
  \item{x}{Realization}
  \item{aic}{The program calls basic R function phi.burg to calculate burg estimates of an AR fit to the data. Aic is turned off and the user specifies the order}
  \item{p}{Order of AR to fit to data}
  \item{ncomp}{Number of additive components to calculate and plot}
}
\value{
  \item{ncomp }{The number of additive components}
  \item{x.comp }{Matrix (i,j) where i designates the component and j denotes time, i.e. (i,j) denotes the ith component at time j}

}
\references{"Applied Time Series Analysis with R, 2nd edition" by Woodward, Elliott, and Gray}
\author{Wayne Woodward}
\examples{data(ss08)
    factor.comp.wge(ss08,p=9,ncomp=4)}
\keyword{ Additive components }
\keyword{ Autoregressive }
