\name{multilm.summaries}
\title{Stabilized Multivariate Tests}
\usage{
summary.multilm(object, test="Hotelling")
}
\alias{summary.multilm}
\arguments{
 \item{object}{an object of class "multilm"}
 \item{test}{the test to be performed: "Hotelling", "SS", "PC-q", "PC-1",
"CS");
}
\description{
 \code{summary.multilm} prints the results of Hotellings T^2 and implements the SS, PC-q, PC-1
and CS test according to Laeuter and Kropf.  
}
\details{
The Hotelling T^2 test is not very useful when the number of observations is
limited but many variables are included in the model. This problem is solved
by the stabilized multivariate test procedures by Laeuter et.al., which are
available in \code{summary.multilm}. All these test procedures follow the same
principle: generate a vector \eqn{d} (or a matrix
\eqn{D}) as a function of \eqn{W = (Y - \bar{Y})^{\top}(Y - \bar{Y}) },
where \eqn{Y} is the matrix of responses. Then analyse the original test problem using
the low dimensional, artificial "responses" \eqn{Yd} or \eqn{YD} by T^2.
The procedures SS, PC-q, PC-1 and CS differ in the way the vector \eqn{d}
or the matrix \eqn{D} is defined. SS, CS and PC-1 test always reduce the
original problem into a one dimensional problem. Naturally, one would use
the univariate procedures instead of a univariate T^2 test. Because of the lack
of a simple standard procedure for testing any linear hypothesis \eqn{H0: K beta =
0}, we use T^2, which reduces to the F-test in the univariate case. 
}
\value{ 
 The test statistic, degrees of freedom and the p-value of the test. 
}
\references{Laeuter, J; Glimm, Ekkehard; Kropf, S. (1998): Multivariate Tests
Based On Left-spherically Distributed Linear Scores, The Annals of
Statistics, Vol. 26, No. 5, pp. 1972-1988}
\author{Torsten Hothorn <hothorn@statistik.uni-dortmund.de>}
\examples{
# Edgar Anderson's famous iris data
data(iris)
# one-classification MANOVA, Y response matrix, X designmatrix
Y <- as.matrix(iris[,1:4]);
x <- c(rep(1,50), rep(0,150), rep(1, 50), rep(0,150), rep(1,50))
X <- matrix(x, ncol=3)
# restrictions: sum of effects equal zero
Z <- c(0,1,1,1);
# test for equal effects
K <- cbind(0,diag(2),-1);
# model (this method returns a multilm object)
mod <- multilm(Y ~ X, K,Z);
# output and stable tests
summary(mod)		# Hotelling T^2: pvalue = 0
summary(mod, "SS")	# SS-Test: pvalue = 0
summary(mod, "PC-q")	# PC-q-Test: pvalue = 0
summary(mod, "PC-1")	# PC-1-Test: pvalue = 0
summary(mod, "CS")	# CS-Test: pvalue = 0
# the iris data is to good -> simulate a "bad" multivariate dataset 
# only 8 observations in 2 groups, 4 variables
observ <- c(1:4,51:54);
rY <- Y[observ,]
rX <- X[observ,1:2]
rZ <- c(0,1,1);
rK <- c(0,1,-1);
rmod <- multilm(rY ~ rX, rK, rZ);
summary(rmod)		# T^2: pvalue = .00052
summary(rmod, "SS")	# SS-Test: pvalue = .008
summary(rmod, "PC-q")	# PC-q-Test: pvalue = .000012
summary(rmod, "PC-1")	# PC-1-Test: pvalue = .0000038
summary(rmod, "CS")	# CS-Test: pvalue = .0002
}
\keyword{multivariate linear models}
