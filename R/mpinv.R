mpinv <- function(X)
{
	# Moore-Penrose Inverse

	# Machine precision

	Eps <- 100 * .Machine$double.eps;
	
	s <- svd(X);
	d <- s$d;
	m <- length(d);
	if (m == 1)
		return(t(s$v%*%(1/d)%*%t(s$u)));

	# remove Eigenvalues equal zero 

	d <- d[d > Eps];
	notnull <- length(d);	
	
	if (notnull == 1)
	{
		inv <- 1/d;
	} else {
		inv <- solve(diag(d));
	}

	# put together again

	if (notnull != m)
	{
		inv <- cbind(inv, matrix(0, nrow=notnull, ncol=(m - notnull)));
		inv <- rbind(inv, matrix(0, nrow=(m-notnull), ncol=m));
	} 
	mp <- s$v%*%inv%*%t(s$u);
 	mp[abs(mp) < Eps] <- 0;
	return(mp);
};
