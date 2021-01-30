\name{plotts.dwt.wge}
\alias{plotts.dwt.wge}
\title{Plots Discrete Wavelet Transform (DWT)}
\description{Plots DWT obtained using functiond dwt from waveslim}
\usage{
plotts.dwt.wge(x, n.levels, type='S8')
}
\arguments{
  \item{x}{Realization (must be of length 2^k for some integer k between 2 and 14}
  \item{n.levels}{Maximum order of discrete wavelet transforms to be calculated. n.levels must be less than or equal to k where n=2^k}
  \item{type}{Discrete wavelet to use: options include 'haar', 'S8','D4','D6',D8'}
}
\details{The wavelsim dwt function names these :'haar', 'la8','d4','d6',and 'd8' respectively and the conversion is done transparently within the R code. This is done transparently within the R code.}
\value{The output is a plot of the DWT.}
\references{Applied Time Series Analysis with R, second edition by Woodward, Gray, and Elliott}
\author{Wayne Woodward}
\note{Requires CRAN package 'waveslim'}
\examples{data(bumps256)
                   plotts.dwt.wge(bumps256,n.levels=4,type='S8')
}
\keyword{  DWT }
\keyword{ discrete wavelet transform }
\keyword{ plot }
