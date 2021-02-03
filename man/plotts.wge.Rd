\name{plotts.wge}
\alias{plotts.wge}
\title{Plot a time series realization}
\description{Given a realization contained in a vector, this function plots it as a time series realization}
\usage{plotts.wge(x,Presentation = 0, xlab = "", ylab = "")}
\arguments{
  \item{x}{The vector containing the time series realization to be plotted}
  \item{Presentation}{If Presentation is 0 then a simple plot of the realization is rendered.  If Presentation is 1 then a ggplot is rendered.}
  \item{xlab}{A string that represents the x-axis label.}
  \item{ylab}{A string that represents the y-axis label.}
}
\value{Simply a plot of the realization}
\references{"Applied Time Series Analysis with R, 2nd edition" by Woodward, Gray, and Elliott}
\author{Wayne Woodward}
\examples{plotts.wge(rnorm(100))}

\keyword{ Plot }
\keyword{ Realization}
