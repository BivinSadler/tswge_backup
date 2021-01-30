\name{hilbert.wge}
\alias{hilbert.wge}
\title{Function to calculate the Hilbert transformation of a given real valued signal(even length)}
\description{Function is used with the tswge function wv.wge }
\usage{
hilbert.wge(input)
}
\arguments{
  \item{input}{realization to be analyzed}
}
\value{
  \item{ans}{Hilbert transformation of the input}
}
\references{Applied Time Series Analysis with R, second edition by Woodward, Gray, and Elliott}
\author{Wayne Woodward}
\examples{
data(airline)
hilbert.wge(airline)
}
\keyword{Hilbert }
\keyword{Wigner-Ville}
