summary.multilm <- function(object, test="Hotelling")
{
	if (!(inherits(object,"multilm"))) stop("no multilm object given");
	z <- .Alias(object)
	Y <- z$response
	n <- nrow(Y)
	X <- z$dmatrix
	K <-z$testK
	Z <- z$restrZ
	ans <- list(method=test, pvalue = 1, stat=0, df1=0, df2=0);
	names(ans$pvalue) <- "P-value";
	if(test=="Hotelling") 
	{ 
		ans$pvalue <- z$hotelp; 
		ans$stat <- z$hotelstat;
		ans$df1 <- z$df1
		ans$df2 <- z$df2
	}
	
	qY <- 1/n*t(Y)%*%rep(1,n)%*%t(rep(1,n))
        W <- (t(Y) - qY)%*%t(t(Y) - qY);
		
	if(test=="PC-q")
	{
	       	# eigenvalue problem: WD = diag(W)DL 
		# diag(W) is non-singular 
		# <=> solve(diag(W)) W D = D L 
		# <=> eigen(solve(diag(W))%*%W)
		# e.g. Numerical Recipes in C, p. 462
                                                       
        	ewp <- eigen(solve(diag(diag(W)))%*%W);
                # Kaiserkriterium
                q <- length(ewp$values[ewp$values > 1]);
                D <- ewp$vectors[,1:q];
		nY <- Y%*%D;
		dummy <- multifit( nY, X, K, Z);
		ans$pvalue <- dummy$hotelp;
		ans$stat <- dummy$hotelstat;
		ans$df1 <- dummy$df1
		ans$df2 <- dummy$df2
	}
	if(test=="PC-1")
	{
		# eigenvalue problem Wd = lambda diag(W) d, see PC-q
                                                          
                ewp <- eigen(solve(diag(diag(W)))%*%W);
                d <- ewp$vectors[,1];
		nY <- Y%*%d
		dummy <- multifit(nY, X, K, Z)
		ans$pvalue <- dummy$hotelp;
		ans$stat <- dummy$hotelstat;
		ans$df1 <- dummy$df1
		ans$df2 <- dummy$df2
	}
	if (test=="SS")
	{
		d <- 1/sqrt(diag(W))
		nY <- Y%*%d
		dummy <- multifit(nY, X, K, Z)
		ans$pvalue <- dummy$hotelp;
		ans$stat <- dummy$hotelstat;
		ans$df1 <- dummy$df1
		ans$df2 <- dummy$df2
	}
	if (test=="CS")
	{
		p <- ncol(W)
		A <- diag(diag(W))
		d <- solve(A)%*%W%*%sqrt(solve(A))%*%rep(1,p) 		
		nY <- Y%*%d
		dummy <- multifit(nY, X, K, Z)
		ans$pvalue <- dummy$hotelp;
		ans$stat <- dummy$hotelstat;
		ans$df1 <- dummy$df1
		ans$df2 <- dummy$df2
	}
	class(ans) <- "summary.multilm" 
	ans
}




