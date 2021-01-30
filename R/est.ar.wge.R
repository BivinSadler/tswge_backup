est.ar.wge=function(x,p=2,factor=TRUE,type="mle"){
method=type
x=x-mean(x)
if (method=="burg") {arest=ar.burg(x,order.max=p,aic=FALSE)
                                    phi.est=arest$ar}
if (method=="yw") {arest=ar.yw(x,order.max=p,aic=FALSE)
                                   phi.est=arest$ar}
if (method=="mle") {arest=est.arma.wge(x,p=p,q=0,factor=TRUE)
                                    phi.est=arest$phi}
if (factor=='TRUE') {factor.wge(phi=phi.est)}
res=backcast.wge(x,phi=phi.est,n.back=50)
avar=0
n=length(x)
for (i in 1:n) {avar=avar+res[i]*res[i]}
avar=avar/n
phi=phi.est
aic=log(avar)+2*(p+1)/n
bic=log(avar)+(p+1)*log(n)/n
aicc=log(avar)+(n+p+1)/(n-p-3)
out1=list(phi=phi,res=res,avar=avar,aic=aic,aicc=aicc,bic=bic)
return(out1)
}

