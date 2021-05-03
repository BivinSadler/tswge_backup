fore.sigplusnoise.wge=function(x,linear=TRUE,method='mle',freq=0,max.p=5,n.ahead=10,lastn=FALSE,plot=TRUE,alpha=.05,limits=TRUE)
{
# if linear=TRUE then a linear trend is fit to the data and the residuals from the trend are fitted usig an AR model
# if linear=FALSE then a cosine function is fit to the data and the residuals from the cosine are fitted usig an AR model
#lastn=TRUE indicates that the last n data values are to be forecast
# lastn=FALSE (default) indicates we want foreacts for n values beyond the end of the realization
n=length(x)
np1=n+1
npn.ahead=n+n.ahead
resid=rep(0,npn.ahead)
xd=rep(0,npn.ahead)
xhat=rep(0,npn.ahead)
pi=3.14159
const=1
xar=rep(0,npn.ahead)
tl=1:n
#
#
if(linear=='TRUE') {
ftl=slr.wge(x)
for(t in 1:n){
xar[t]=x[t]-ftl$b0hat-t*ftl$b1hat
xd[t]=x[t]-xar[t]
                   }
#
if(lastn=='FALSE') {for(t in np1:npn.ahead) {xd[t]=ftl$b0hat+t*ftl$b1hat}}
}
#
if(linear=='FALSE') {
x1=rep(0,n)
x2=rep(0,n)
for(t in 1:n) {
x1[t]=cos(2*pi*freq*t)
x2[t]=sin(2*pi*freq*t)}
xm=rbind(x1,x2)
xmt=t(xm)
ftc=lm(x~xmt)
#
for(t in 1:n){
xar[t]=x[t]-ftc$coefficients[1]-ftc$coefficients[2]*x1[t]-ftc$coefficients[3]*x2[t]
xd[t]=x[t]-xar[t]
                   }
if(lastn=='FALSE') {
x1=rep(0,npn.ahead)
x2=rep(0,npn.ahead)
for(t in 1:npn.ahead){
x1[t]=cos(2*pi*freq*t)
x2[t]=sin(2*pi*freq*t)
xd[t]=ftc$coefficients[1]+ftc$coefficients[2]*x1[t]+ftc$coefficients[3]*x2[t]
}
}
}
order=aic.wge(xar[1:n],p=0:max.p,q=0:0)
p=order$p
phi=0
method.est=method
if(p > 0) {w=est.ar.wge(xar[1:n],p=p,method=method.est)
phi=w$phi
}
#
if (p > 0) {for(jp in 1:p) {const=const-phi[jp]}}
#
#
# Calculate Box-Jenkins Forecasts
#
#
#Calculating Residuals
#
#
resid=backcast.wge(xar[1:n],phi,theta=0,n.back=50)

p1=p+1
xbar=mean(xar)
maconst=const*xbar
#for (i in p1:n) {resid[i]=xar[i]
#   if ( p > 0) {for (jp in 1:p) {resid[i]=resid[i]-phi[jp]*xar[i-jp]}
#                   resid[i]=resid[i]-maconst}}
#
# Calculating Forecasts for AR noise
#
#
npn.ahead=n+n.ahead
xhat=rep(0,npn.ahead)
mm=n
#
if(lastn==TRUE) {mm=n-n.ahead}
#
for (i in 1:mm) {xhat[i]=xar[i]}
for (h in 1:n.ahead) {
if (p > 0) {for (jp in 1:p) {xhat[mm+h]=xhat[mm+h]+phi[jp]*xhat[mm+h-jp]}
                    xhat[mm+h]=xhat[mm+h]+maconst}
                                 }
# xhat's are forecasts for noise - they need to have the signal added back in
#
#   Calculate psi weights for forecasts limits for AR noise forecasts limits
#
#
xi=psi.weights.wge(phi,theta=0,lag.max=n.ahead)
#
#
#
#Setting up for plots
nap1=n.ahead+1
fplot=rep(0,nap1)
maxh=mm+n.ahead
llplot=rep(0,nap1)
ulplot=rep(0,nap1)
f=rep(0,nap1)
ll=rep(0,nap1)
ul=rep(0,nap1)
wnv=0
xisq=rep(0,n.ahead)
se=rep(0,n.ahead)
se0=1
for (i in 1:n) {wnv=wnv+resid[i]**2}
wnv=wnv/n
xisq[1]=1
for (i in 2:n.ahead) {xisq[i]=xisq[i-1]+xi[i-1]^2}
for (i in 1:n.ahead) {se[i]=sqrt(wnv*xisq[i])}
#
#
# calculate forecasts for signal+noise
#
#

fplot[1]=x[mm]
for (i in 1:n.ahead) {fplot[i+1]=xhat[mm+i]+xd[mm+i]}
ulplot[1]=x[mm]
for (i in 1:n.ahead) { ulplot[i+1]=fplot[i+1]+1.96*se[i]}
llplot[1]=x[mm]
for (i in 1:n.ahead) { llplot[i+1]=fplot[i+1]-1.96*se[i]}
#
if(limits==FALSE) {
  if(lastn==TRUE) {maxp=max(x,xhat[1:n]+xd[1:n])
  minp=min(x,xhat[1:n]+xd[1:n])
}
  else            {maxp=max(x,fplot)
  minp=min(x,fplot)}}


if(limits==TRUE) {
maxp=max(x,ulplot)
 minp=min(x,llplot)
}
numrows <- 1
numcols <- 1
timelab <- 'Time'
valuelab <- ''
fig.width <- 5
fig.height <- 2.5
cex.labs <- c(.8,.7,.8)
par(mfrow=c(numrows,numcols),mar=c(3.8,2.5,1,1))
t<-1:n;
np1=n+1
np.ahead=mm+n.ahead
tf<-mm:np.ahead
if (plot=='TRUE') {
fig.width <- 5
fig.height <- 2.5
cex.labs <- c(.8,.7,.8)
par(mfrow=c(numrows,numcols),mar=c(3.8,2.5,1,1))
plot(t,x,type='o',xaxt='n',yaxt='n',cex=.8,pch=16,cex.lab=.75,cex.axis=.75,lwd=.75,xlab='',ylab='',xlim=c(1,maxh),ylim=c(minp,maxp))
axis(side=1,cex.axis=.8,mgp=c(3,0.15,0),tcl=-.3);
axis(side=2,las=1,cex.axis=.8,mgp=c(3,.4,0),tcl=-.3)
mtext(side=c(1,2,1),cex=cex.labs,text=c(timelab,valuelab,""),line=c(.8,1.1,1.8))
points(tf,fplot,type='o',lty=3,cex=1,lwd=2,pch=18);
#points(t,xd)
if(limits=='TRUE') {points(tf,ulplot,type='o',lty=3,cex=0.6,lwd=.75,pch=18)
points(tf,llplot,type='o',lty=3,cex=0.6,lwd=.75,pch=18) }
    }
np1=n+1
nap1=n.ahead+1
f=fplot[2:nap1]
ll=llplot[2:nap1]
ul=ulplot[2:nap1]
if (linear==TRUE){
out1=list(b0hat=ftl$b0hat,b1hat=ftl$b1hat,sig=xd,z=xar[1:n],phi.z=phi,f=f,ll=ll,ul=ul,resid=resid,wnv=wnv,se=se,xi=xi)
return(out1)}
if (linear==FALSE){
b0=ftc$coefficients[1]
b1=ftc$coefficients[1]
b2=ftc$coefficients[3]
out1=list(b=ftc$coefficients,sig=xd,z=xar[1:n],phi.z=phi,f=f,ll=ll,ul=ul,resid=resid,wnv=wnv,se=se,xi=xi)
return(out1)}
}


