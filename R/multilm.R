multilm <- function(formula, K=0 , Z=0, data=list())
{
	mf <- model.frame(formula, data=data)
        Y <- model.response(mf);
        X <- model.matrix(formula, data=data);
	return(multifit(Y,X,K,Z));
}
		 
	
	