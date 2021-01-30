\name{mult.wge}
\alias{mult.wge}
\title{Multiply Factors}
\description{The function multiplies the AR (or MA) factors of a model to produce the model in unfactored form. Requires the CRAN package 'PolynomF'.}
\usage{mult.wge(fac1 = 0, fac2 = 0, fac3 = 0, fac4 = 0, fac5 = 0, fac6 = 0)
}
\arguments{
  \item{fac1}{First factor to be multiplied}
  \item{fac2}{Second factor to be multiplied}
  \item{fac3}{Third factor to be multiplied (you may use a maximum of 6 factors)}
  \item{fac4}{Fourth factor to be multiplied (you may use a maximum of 6 factors)}
  \item{fac5}{Fifth factor to be multiplied (you may use a maximum of 6 factors)}
  \item{fac6}{Sixth factor to be multiplied (you may use a maximum of 6 factors)}
}
\value{
  \item{char.poly }{The characteristics polynomial of the full model}
  \item{model.coef }{Model coefficients of the full model using notation in "Applied Time Series Analysis, 2nd edition" by Woodward, Gray, and Elliott}
}
\references{"Applied Time Series Analysis with R, 2nd edition" by Woodward, Gray, and Elliott}
\author{Wayne Woodward}
\note{Requires CRAN package 'PolynomF'}
\examples{
fac1=c(1.6,-.9)
fac2=.8
mult.wge(fac1,fac2)
}
\keyword{ Factors }

