est.arma.wge=function (x, p = 0, q = 0,factor=TRUE) 
{
    xbar=mean(x)
    x = x - mean(x)
    w = getOption("warn")
    options(warn = -1)
    a = arima(x, c(p, 0, q))
    options(warn = w)
    c = as.vector(coef(a))
    v = as.vector(diag(vcov(a)))
    if (p == 0) {
        phi = 0
        se.phi = 0
    }
    else {
        phi = c[1:p]
        se.phi = v[1:p]
    }
    if (q == 0) {
        theta = 0
        se.theta = 0
    }
    else {
        theta = -c[(p + 1):(p + q)]
        se.theta = v[(p + 1):(p + q)]
    }

if (factor==TRUE & p>0) {factor.wge(phi=phi,theta=theta)}
res=backcast.wge(x,phi=phi,theta=theta,n.back=50)
avar=0
n=length(x)
for (i in 1:n) {avar=avar+res[i]*res[i]}
avar=avar/n
aic=log(avar)+2*(p+q+1)/n
bic=log(avar)+(p+q+1)*log(n)/n
aicc=log(avar)+(n+p+q+1)/(n-p-q-3)
out1 = list(phi = phi, theta = theta, res=res,avar=avar, xbar=xbar,aic=aic,aicc=aicc,bic=bic, 
        se.phi = se.phi, se.theta = se.theta)
#    a = .innovation.update(x, a)
    return(out1)
}
