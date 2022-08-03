\name{unit.circle.wge}
\alias{unit.circle.wge}
\title{Plot the roots of the characteristic equation on the complex plain.}
\description{This function plots the roots of the characteristic equation on the complex plain and super imposes the Unit Circle to show if a root is inside, outside or on the Unit Circle.  The modulus and absolule reciprical are also displayed. }
\usage{unit.circle.wge(real = 0, imaginary = 0)
}
\arguments{
  \item{real}{the real part of the root}
  \item{imaginary}{the imaginary part of the root}
}
\value{returns a plot of the root with respect to the unit circle}
\references{"Applied Time Series Analysis with R, 2nd edition" by Woodward, Gray, and Elliott}
\author{Bivin Sadler}
\examples{unit.circle.wge = function(real = .9, imaginary = .95)}
\keyword{ Unit Circle }
\keyword{ Root}
\keyword{ Stationary}
