\name{bitcoin}
\alias{bitcoin}
\docType{data}
\title{Daily Bitcoin Prices From Jan 1, 2020 to April 5, 2021}
\description{This dataset contains the daily price of bitcoin from January 1, 2021 to April 5, 2021. The data was gather from Yahoo Finance and included missing values on April 17 and October 9, 12 and 13 of 2020.  These values were imputed using the mean of the values immediatly before and after them (linear interpolation).}
\usage{data("bitcoin")}
\format{
  The format is:
 num [1:461] 7200.174 6985.470 7344.884  ...
}
\source{Yahoo Finance}
\references{"The Joy of Forecasting" by Woodward, Gray, and Elliott}
\examples{data(bitcoin)}
\keyword{datasets}
