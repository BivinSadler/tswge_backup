gen.arch.wge=function(n,alpha0,alpha, plot='TRUE',sn=0)
{
# if sn=0 then produces random results, otherwise, uses seed >0
  if (sn > 0) {set.seed(sn)}
q=length(alpha)
# setting spin length to 1000
spin=1000
ntot=n+spin+q
eps<-rnorm(ntot, mean=0, sd=1);
sig2 <- rep(1, ntot)
asq=rep(0,ntot)
qp1=q+1
for (t in q) { sig2[t]=1 }
for(t in qp1:ntot) {
      sig2[t]<-alpha0
       for (i in 1:q) {sig2[t]=sig2[t]+alpha[i]*asq[t-i]^2}
       asq[t]=eps[t]*(sig2[t])^.5}
spinq=spin+q
spinq1=spinq+1
real=asq[spinq1:ntot]
if(plot== 'TRUE') {plot(real,type='l',xlab='')}
return(real)
}


