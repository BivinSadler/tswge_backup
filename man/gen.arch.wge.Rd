\name{gen.arch.wge}
\alias{gen.arch.wge}
\title{Generate a realization from an ARCH(q0) model}
\description{Generates a realization of length n from the GARCH(q0) model (4.23) in "Applied Time Series Analysis with R, 2nd edition" by Woodward, Gray, and Elliott}
\usage{
gen.arch.wge(n, alpha0, alpha, plot = TRUE,sn=0)
}
\arguments{
  \item{n}{Length of realization to be generated}
  \item{alpha0}{The constant alpha0 in model (4.23)}
  \item{alpha}{A vector of length q0 containing alpha1 through alphaq0}
  \item{plot}{If plot=TRUE (default) the generated realization is plotted}
  \item{sn}{determines the seed used in the simulation. sn=0 produces new/random realization each time. sn=positive integer produces same realization each time}
}
\value{returns the generated realization}
\references{"Applied Time Series Analysis with R, 2nd edition" by Woodward, Gray, and Elliott}
\author{Wayne Woodward}
\examples{gen.arch.wge(n=200,alpha0=.1,alpha=c(.36,.27,.18,.09))}
\keyword{ ARCH }
\keyword{ Conditional variance}
