\name{mpinv}
\title{Moore-Penrose-Inverse}
\usage{
mpinv(X)
}
\alias{mpinv}
\arguments{
 \item{X}{a matrix, which Moore-Penrose-Inverse is to calculate}
}
\description{
 \code{mpinv} returns the Moore-Penrose-Inverse. 
}
\references{any Linare Algebra book}
\author{Torsten Hothorn <hothorn@statistik.uni-dortmund.de>}
\examples{
X <- cbind(1, diag(10))	# reduced rank
mpinv(X) 
}
\keyword{linear models}