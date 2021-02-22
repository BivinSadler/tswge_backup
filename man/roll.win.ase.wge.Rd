\name{roll.win.wge}
\alias{roll.win.wge}
\title{Function to Calculate the Rolling Window ASE for Variable Window Sizes}
\description{This function creates as many "windows" as is possible wiht the data and calculates an ASE for each window.  The resulting "rolling window ASE" is the average of the individual ASEs from each window. Specfically, a "window" is a subset of the data of size (traininSize + horizon)."}
\usage{
roll.win.ase.wge = function(series, trainingSize, horizon = 1, s = 0, d = 0, phis = 0, thetas = 0)
}
\arguments{
  \item{series}{The data}
  \item{horizon}{The number of observations ahead to be forecasted.}
  \item{s}{Order of the seasonal difference, default=1}
  \item{d}{Order of the difference}
  \item{phis}{Vector of AR coefficients}
  \item{thetas}{Vector of MA coefficients}
}
\value{
  \item{rwASE}{The average of the individual ASEs of each window}
  \item{numwindows}{The number of windows}
  \item{horizon}{The number of observations ahead to be forecasted.}
  \item{s}{Order of the seasonal difference, default=1}
  \item{d}{Order of the difference}
  \item{phis}{Vector of AR coefficients}
  \item{thetas}{Vector of MA coefficients}
}
\references{"The new book."}
\author{Bivin Sadler}

\keyword{ rolling }
\keyword{ window }
\keyword{ ASE }