print.summary.multilm <- function (x) 
{
    cat("test procedure: ", x$method, "\n")
    cat("test statistic: ", x$stat, "\n")
    cat("degrees of freedom DF1: ", x$df1, " DF2: ", x$df2, "\n")
    cat("p-value: ", x$pvalue, "\n")
    invisible(x)
}
