\name{wv.wge}
\alias{wv.wge}
\title{Function to calculate Wigner Ville spectrum}
\description{Calculates and plots Wigner-Ville spectrum for a realization}
\usage{
wv.wge(x)
}
\arguments{
  \item{x}{Realization to be analyzed}
}
\value{Plots Wigner-Ville spectrum}
\references{Boashash (2003). Time Frequency Analysis}
\author{Wayne Woodward}
\examples{
data(doppler)
wv.dop=wv.wge(doppler)
}
\keyword{ Wigner-Ville spectrum }
\keyword{ TVF}
