\name{roll.win.wge}
\alias{roll.win.wge}
\title{Function to Calculate the Rolling Window ASE for Variable Window Sizes}
\description{This function creates as many "windows" as is possible wiht the data and calculates an ASE for each window.  The resulting "rolling window ASE" is the average of the individual ASEs from each window. Specfically, a "window" is a subset of the data of size (traininSize + horizon)."}
\usage{
roll.win.ase.wge = function(series, trainingSize, horizon = 1, s = 0, d = 0, phis = 0, thetas = 0)
}
\arguments{
  \item{series}{The data}
  \item{d}{Order of the difference}
  \item{s}{Order of the seasonal difference, default=1}
  \item{trainingSize}{Size of the data used to generage the forecast.}
  \item{horizon}{The number of observations ahead to be forecasted.}
  \item{phi}{Vector of AR coefficients}
  \item{theta}{Vector of MA coefficients}
}
\value
  \item{rwASE} the average of the individual ASEs of each window.
  
\references{"The new book."}
\author{Bivin Sadler}

\keyword{ rolling }
\keyword{ window }
\keyword{ ASE }