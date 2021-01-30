
#
#
factor.comp.wge=function(x,aic=FALSE,p,ncomp,report="c:\\ATSA Reports\\report.txt",append=TRUE)
{
#
#x is the time series realization to be analyzed#
# order is order of the AR model to be fit to data if aic=FALSE
#       is max AR order if aic=TRUE
#
#
#cat('Coefficients of Original polynomial: ', phi,'\n')
#cat(file=report,'Coefficients of Original polynomial: ', phi,'\n')
#cat(file=report,append=app,'\n')
n=length(x)
phi.mle=ar.mle(x,order.max=p,aic=FALSE)
phi=phi.mle$ar
phi
factor.wge(phi,report)
mphi=-phi
one=c(1)
coef=c(one,mphi)   # these are the actual coefficients of the pth order polynomial
root1=polyroot(coef)    # with signs in characteristic equation
roota=abs(root1)
perm=sort(roota,index.return=TRUE)
nr=length(root1)
title=rep(0,ncomp)
roots=rep(0,nr)
for(i in 1:nr) {ii=perm$ix[i]
                roots[i]=root1[ii]
               }
root1=roots
nfactors <- 0.5*(length(root1)+sum(abs(Im(root1)) <= 10^(-5)))
if (ncomp > nfactors) {ncomp=nfactors}
nf1=nfactors+1
#par(mfrow=c(5,1),mar=c(3.8,2.5,1,1))#
maxjcomp=rep(0,nf1)
fac=matrix(0,nr,nr)
for (i in 1:nr) {fac[i,]=(1/root1)^(nr-i)}
fac
facinv=solve(fac)
coefs.fac=matrix(0,nfactors,p)
#***get nfactors
#  k=# components (<=nfactors)
j=0
for (i in 1:nfactors) {if(abs(Im(root1[j+1])) <= 10^(-5)) 
                {coefs.fac[i,]=facinv[j+1,]
                 jump=1
#          
                }
                else {coefs.fac[i,]=facinv[j+1,]+facinv[j+2,]
                 jump=2
                                }
j=j+jump}

#for (j in 1:ncomp) {cat(j,coefs.fac[j,],'\n') }
#cat(coefs.fac[1,],coefs.fac[2,],'\n')
x.comp=matrix(0,ncomp,n)
#
for (jcomp in 1:ncomp) {
for (t in p:n) {
for (j in 1:p) {x.comp[jcomp,t]=x.comp[jcomp,t]+coefs.fac[jcomp,j]*x[t-p+j]}
                }
                           }
x.comp=Re(x.comp)
#for (j in 1:nfactors) {maxjcomp(jcomp)=max(x.comp[j,])}
#for (j in 1:nfactors) {minjcomp(jcomp)=min(x.comp[j,])}                         
maxcomp=max(x.comp)
mincomp=min(x.comp)
maxx=max(x)
minx=min(x)
max.tot=max(maxx,maxcomp)
min.tot=min(minx,mincomp)

ncomp1=ncomp+1
par(mfrow=c(ncomp1,1),mar=c(3.8,2.5,1,2))
t=1:n
plot(t,x,type='l',ylim=c(min.tot,max.tot),xaxt='n',yaxt='n',cex=0.8,cex.lab=.75,cex.axis=.75,lwd=.75,xlab='',ylab='',pch=16)
axis(side=1,cex.axis=1,mgp=c(3,0.15,0),tcl=-.3);
axis(side=2,las=1,cex.axis=1,mgp=c(3,.4,0),tcl=-.3)
mtext(side=c(1,2,1),cex=cex.labs,text=c('Time','','Realization'),line=c(1.1,1.3,2),font=c(1,1,1))
t=p:n
for (j in 1:ncomp) {x.c=x.comp[j,]
#cat('j',j,'\n')
#cat(x.c,'\n')
titl=sprintf('Component %d',j)
plot(t,x.c[p:n],type='l',ylim=c(min.tot,max.tot),xaxt='n',yaxt='n',cex=0.8,cex.lab=.75,cex.axis=.75,lwd=.75,xlab='',ylab='',pch=16)
axis(side=1,cex.axis=1,mgp=c(3,0.15,0),tcl=-.3);
axis(side=2,las=1,cex.axis=1,mgp=c(3,.4,0),tcl=-.3)
mtext(side=c(1,2,1),cex=cex.labs,text=c('Time','',titl),line=c(1.1,1.3,1.8),font=c(1,1,1))
                       }

out1=list(ncomp=ncomp,x.comp=x.comp)
return(out1)   
}
#
#