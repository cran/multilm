\name{summary.multilm}
\title{Stable Multivariate Tests}
\usage{
summary.multilm(object, test="Hotelling")
}
\alias{summary.multilm}
\arguments{
 \item{object}{an object of class "multilm"}
 \item{test} {the test to be performed: "Hotelling", "SS", "PC-q", "PC-1",
"CS");
}
\description{
 \code{summary.multilm} implements the PC-q and PC-1 test according to
Laeuter and Kropf (1996)  
}
\details{
The Hotelling \eqn{T^2} test is not very useful when the number of observations is
limited but many variables are included in the model. This problem is solved
by the stable multivariate procedures by Laeuter and Kropf (1996), which are
available in \code{summary.multilm}. All these test procedures follow the same
principle: generate a score vectors \code{d} (or a matrix of score vectors
\code{D}) as a function of \code{W = (Y - \bar{Y})^{\top}(Y - \bar{Y})},
where \code{Y} is the matrix of responses. Then analyse the original test problem using
the lowdimensional, artificial "responses" \code{Yd} or \code{YD} by \eqn{T^2}.
The procedures "SS", "PC-q", "PC-1" and "CS" differ in the way the score
vector or matrix is defined. "SS", "CS" and "PC-1" test always reduce the
original problem into a one dimensional problem. Naturally, one would use
the univariate procedures instead of a univariate \eqn{T^2}. Because of the lack
of a simple standard procedure testing linear hypotheses \eqn{H0: K\beta = 0}, we
use \eqn{T^2}, which reduces to the F-test in the univariate case. 
}
\value{ 
 The p-value of the test. 
}
\references{Laeuter, J; Glimm, Ekkehard; Kropf, S. (1998): Multivariate Tests
Based On Left-spherically Distributed Linear Scores, The Annals of
Statistics, Vol. 26, No. 5, pp. 1972-1988}
\author{Torsten Hothorn <hothorn@amadeus.statistik.uni-dortmund.de>}
\examples{
# Edgar Anderson's famous iris data
data(iris)
# one-classification MANOVA, Y response matrix, X designmatrix
Y <- as.matrix(iris[,1:4]);
X <- cbind(c(rep(1,50), rep(0,100)), c(rep(0,50), rep(1,50), rep(0,50)), c(rep(0,100), rep(1,50)));
# restrictions: sum of effects equal zero
Z <- c(0,1,1,1);
# test for equal effects
K <- cbind(0,diag(2),-1);
# model (this method returns a multilm object)
mod <- multilm(Y ~ X, K,Z);
# output and stable tests
summary.multilm(mod)		# Hotelling T^2: pvalue = 0
summary.multilm(mod, "SS")	# SS-Test: pvalue = 0
summary.multilm(mod, "PC-q")	# PC-q-Test: pvalue = 0
summary.multilm(mod, "PC-1")	# PC-1-Test: pvalue = 0
summary.multilm(mod, "CS")	# CS-Test: pvalue = 0
# the iris data is to good -> simulate a "bad" multivariate dataset 
# only 8 observations in 2 groups, 4 variables
observ <- c(1:4,51:54);
rY <- Y[observ,]
rX <- X[observ,1:2]
rZ <- c(0,1,1);
rK <- c(0,1,-1);
rmod <- multilm(rY ~ rX, rK, rZ);
summary.multilm(rmod)$pvalue		# T^2: pvalue = .00052
summary.multilm(rmod, "SS")	# SS-Test: pvalue = .008
summary.multilm(rmod, "PC-q")	# PC-q-Test: pvalue = .000012
summary.multilm(rmod, "PC-1")	# PC-1-Test: pvalue = .0000038
summary.multilm(rmod, "CS")	# CS-Test: pvalue = .0002
}
\keyword{multivariate linear models}