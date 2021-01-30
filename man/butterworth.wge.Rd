\name{butterworth.wge}
\alias{butterworth.wge}
\title{Perform Butterworth Filter}
\description{The user can specify the order of the filter, and whether it is low pass ("low"), high pass ("high"), band stop ("stop"), or band pass ("pass") filter. Requires the CRAN package 'signal'.}
\usage{butterworth.wge(x, order, type, cutoff,plot=TRUE)}
\arguments{
  \item{x}{Realization to be filtered}
  \item{order}{Order of the Butterworth filter}
  \item{type}{Either "low", "high", "stop", or "pass" as dicsussed in Descriptions}
  \item{cutoff}{For "low" and "high": cutoff is a real number. For "stop" and "band": cutoff is a 2-component vector}
  \item{plot}{If plot=TRUE then plots of the original and filtered data are produced.}
}
\value{The filtered data}
\references{"Applied Time Series Analysis with R, 2nd edition" by Woodward, Gray, and Elliott}
\author{Wayne Woodward}
\note{Requires CRAN package 'signal'}
\examples{data(wages)
    butterworth.wge(wages,order=4,type="low",cutoff=.05)}
\keyword{ filtering }
\keyword{ Butterworth }
