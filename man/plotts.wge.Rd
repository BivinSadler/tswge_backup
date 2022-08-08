\name{plotts.wge}
\alias{plotts.wge}
\title{Plot a time series realization}
\description{Given a realization contained in a vector, this function plots it as a time series realization}
\usage{plotts.wge(x,style = 0, xlab = "Time", ylab = "",main="",col='black',text_size=12,
lwd=0.75,cex=0.5,cex.lab=0.75,cex.axis=0.75,xlim=NULL,ylim=NULL)}
\arguments{
  \item{x}{The vector containing the time series realization to be plotted}
  \item{style}{If style is 0 then a simple plot of the realization is rendered.  If style is 1 then a ggplot is rendered.}
  \item{xlab}{A string that represents the x-axis label.}
  \item{ylab}{A string that represents the y-axis label.}
  \item{main}{A string that represents the main title.}
  \item{col}{Color of plot.}
  \item{text_size}{Text size.}
  \item{lwd}{Line width.}
  \item{cex}{See R documentation.}
  \item{cex.lab}{See R documentation.}
  \item{cex.axis}{See R documentation.}
  \item{xlim}{String giving x-axis plot limits.}
  \item{ylim}{String giving y-axis plot limits.}
}
\value{Simply a plot of the realization}
\references{"Applied Time Series Analysis with R, 2nd edition" by Woodward, Gray, and Elliott}
\author{Wayne Woodward}
\examples{data(sunspot2.0);plotts.wge(sunspot2.0)}
\keyword{ Plot }
\keyword{ Realization}