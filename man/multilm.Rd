\name{multilm}
\title{Mulivariate Linear Models}
\usage{
multilm(formula, K, Z, data=list())
}
\alias{multilm}
\arguments{
 \item{formula}{a symbolic description for the multivariate model to be tested}
 \item{K}{a test matrix}
 \item{Z}{a matrix for parameter restrictions} 
\item{data}{an optional data frame containing the variables in the model.
By default the variables are taken from the environment which 'multilm' is
called from}
}
\description{
 \code{multilm} fits a multivariate linear model und performs the Hotelling
T^2 - Test for a given linear test problem. 
}
\details{
A multivariate linear model is a model of the form \eqn{Y = X B + E}, where
\eqn{Y} is
the matrix of responses, \eqn{X} is the design matrix, \eqn{B} is the matrix of
coefficients and \eqn{E} a matrix of normally distributed errors. Parameter restrictions
can be included by the \eqn{Z} matrix: \eqn{Z B = 0} (which has applications in
MANOVA). \code{multilm} additionally calculates the Hotelling T^2-Test for
the given test problem: \eqn{H0: K B = 0}. An approximation by Laeuter
is used for the distribution of the T^2-statistic (and therefore 
for the p-value). T^2 is not very useful when the number of observations is
limited but many variables are included in the model. This problem is solved
by the stabilized multivariate test procedures by Laeuter et. al. , which are
available in \code{summary.multilm}  
}
\value{
 A object of class \code{multilm} containing the following components:
  \item{hotelstat}{the T^2 test statistic}
  \item{hotelp}{the pvalue of the T^2-test}
  \item{coefficients}{the matrix of estimated coefficients}
  \item{covar}{the estimation of the covariance matrix} 
}
\references{Ahrens, H.; Laeuter, J. (1981): Mehrdimensionale Varianzanalyse,
Berlin}
\author{Torsten Hothorn <hothorn@statistik.uni-dortmund.de>}

\seealso{\code{\link{summary.multilm}}}

\examples{
# Edgar Anderson's famous iris data
data(iris)
# one-classification MANOVA, Y response matrix, X design matrix
Y <- as.matrix(iris[,1:4]);
x <- c(rep(1,50), rep(0,150), rep(1, 50), rep(0, 150), rep(1,50))
X <- matrix(x, ncol=3)
# restrictions: sum of effects equal zero
Z <- c(0,1,1,1);
# test for equal effects
K <- cbind(0,diag(2),-1);
# model (this method returns a multilm object)
mod <- multilm(Y ~ X, K,Z);
# output and stable tests
summary(mod)		# Hotelling T^2: pvalue = 0
}
\keyword{multivariate linear models}
