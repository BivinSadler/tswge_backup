gen.garch.wge=function(n,alpha0,alpha,beta,plot=TRUE,sn=0)
{
#n is the length of the realization to be generated
#alpha0 is the value alpha0 given in (4.26)
#alpha is a vector of length q0 containing the alpha coefficients
#beta is a vector of length p0 containing the beta coefficients
# setting spin length to 1000
q0=length(alpha)
p0=length(beta)
# if sn=0 then produces random results, otherwise, uses seed >0
  if (sn > 0) {set.seed(sn)}
spin=1000
ntot=n+spin+q0+p0
#
eps<-rnorm(ntot, mean=0, sd=1);
sig2 <- rep(1, ntot)
asq=rep(0,ntot)
qp1=q0+1
for (t in q0) { sig2[t]=1 }
for(t in qp1:ntot) {
      sig2[t]<-alpha0
       for (i in 1:q0) {sig2[t]=sig2[t]+alpha[i]*asq[t-i]^2}
       for (j in 1:p0) {sig2[t]=sig2[t]+beta[j]*sig2[t-j]}
       asq[t]=eps[t]*(sig2[t])^.5}
spinq=spin+q0+p0
spinq1=spinq+1
real=asq[spinq1:ntot]
if(plot== 'TRUE') {plot(real,type='l')}
return(real)
}


