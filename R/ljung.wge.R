ljung.wge=function (x, K=24, p=0,q=0) 
{
    cor <- acf(x, lag.max = K, plot = FALSE, na.action = na.pass)
    n <- length(x)
    df=K-p-q
    obs <- cor$acf[2:(K + 1)]
cat('Obs',obs,'\n')
    test <- "Ljung-Box test"
    chi.square <- n * (n + 2) * sum(1/seq.int(n - 1, n - K) * obs^2)
    pval <- 1 - pchisq(chi.square, K - p-q)
out1=list(test=test,K=K, chi.square=chi.square,df=df,pval=pval)
return(out1)   
}
