multifit <- function(Y,X,K,Z)
{

	if (is.vector(Y)) 
	{ 
		# multivariate T^2 reduces to F-Test in the univariate case
		n <- length(Y)
		p <- 1
		Y <- matrix(Y, ncol=1)
	}

	z <- list(response = Y, dmatrix = X, testK = K, restrZ = Z, coefficients = 0, covar = 0, df1 = 0, df2 = 0, hotelstat = 0 , hotelp = 1); 
	class(z) <- "multilm"

	# multivariate linear model Y = X Phi + E

	n <- nrow(Y);
	p <- ncol(Y);
	m <- ncol(X);
	
	# if no test matrix specified, test: all coeffcients equal 0

	if (all(K == 0)) K <- cbind(diag(m-1), -1)
	
	# Kern(Z)-rest. Moore-Penrose-Inverse of X

	MPXZ <- rmpinv(X,Z);
	
	if (all(Z == 0))
	{
		ZK <- K;
	} else {
		ZK <- rbind(Z,K);
	}
	MPXZK <- rmpinv(X,ZK)
	
	# parameter estimation

	hatPHI <- MPXZ%*%Y;
	
	z$coefficients <- hatPHI;

	# projection matrices

	Ph <- X%*%(MPXZ - MPXZK);
	Pe <- diag(n) - X%*%MPXZ;
	Sh <- t(Y)%*%Ph%*%Y;
	Se <- t(Y)%*%Pe%*%Y;
	
	# degrees of freedom

	nh <- sum(diag(Ph));
	ne <- sum(diag(Pe));
	
	# covariance estimation 

	hatSigma <- 1/ne*Se;

	z$covar <- hatSigma

	# Hotelling-Lawley Test: H0: K%*%Phi = 0
	
	eigenvalues <- eigen(Sh%*%solve(Se))$val  

	realeigenw <- Re(eigenvalues[Im(eigenvalues) == 0]);

	HL <- sum(realeigenw);

	z$hotelstat <- HL

	# approximation according to Laeuter (1974): Approximation des
	# Hotellingschen T^2 durch die F-Verteilung, Biometrische 
	# Zeitschrift 16, Heft 3 
	
	f1 <- nrow(K);
	if (is.null(f1)) f1 <- 1;
	f2 <- nrow(Y) - (ncol(X) - 1);		# Rang(X) = m - 1
	p <- ncol(Y);
	if (f1 + f2 - f1*p -1 > 0)
	{
		g1 <- (f1*p*(f2 - p))/(f1 + f2 - f1*p - 1)
		g2 <- f2 - p +1
	} else {
		z$hotelp <- 1
		stop("cannot approximate null distribution: to less observations");
	}
	tildeF <- (f2 - p -1)/(f1*p)*g2/(g2 -2)*HL;
	z$hotelstat <- tildeF
	z$df1 <- round(g1)
	z$df2 <- round(g2)
	pvalue <- 1 - pf(tildeF, round(g1), round(g2));

	z$hotelp <- pvalue

	return(z)
}

	
		 
	
	