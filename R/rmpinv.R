rmpinv <- function(X,Z = 0)
{
	# Kern(Z)-restricted Moore-Penrose-Inverse
	
	# Machinen precision

	Eps <- 100* .Machine$double.eps;
	
	
	if (all(Z == 0))
	{
		return(mpinv(X));
	} else {
		ZZ <- mpinv(Z)%*%Z;
		ZZ[abs(ZZ) < Eps] <- 0;
		n <- nrow(ZZ);
		if (all(diag(n) == ZZ))
		{
			stop("Error: ZZ equals diag(n)");
			return(-1);
		} else {
			DZ <- diag(n) - ZZ;
			DZ[abs(DZ) < Eps] <- 0;
			return(mpinv(X%*%DZ));
		}
	}
}
