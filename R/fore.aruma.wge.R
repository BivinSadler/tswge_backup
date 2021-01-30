fore.aruma.wge=function(x,phi=0,theta=0,d=0,s=0,lambda=0, n.ahead=2,lastn=FALSE, plot=TRUE,limits=TRUE)
{
n=length(x)
p=length(phi)
if(sum(phi^2)==0) {p=0
            fac1=0}
else if(sum(phi^2)>0) {fac1=phi}
q=length(theta)
if(sum(theta^2)==0) {q=0}
resid=rep(0,n)
npn.ahead=n+n.ahead
xhat=rep(0,npn.ahead)
xbar=mean(x)
if (d > 3) {cat('d>3','\n')
            return(d)}
if (d==0) {diffac=0}
else if (d==1) {diffac=1}
else if (d==2) {diffac=c(2,-1)}
else if (d==3) {diffac=c(3,-3,1)}
if (d<4) {fac2=diffac} 
#cat('s=',s,'\n') 
if(s==0) {fac3=0}
if(s > 0) {seas=rep(0,s)
           seas[s]=1
           fac3=seas}
#cat('seas',seas,'\n')
dlam=length(lambda)
if(sum(lambda^2)==0) {dlam=0
            fac4=0}
else if(sum(lambda^2)>0) {fac4=lambda}
#cat('lambda',lambda,'\n')
prod.tmp=mult.wge(fac1,fac2,fac3,fac4)
phitot=prod.tmp$model.coef
#cat('fac1',fac1,'\n')
#cat('fac2',fac2,'\n')
#cat('fac3',fac3,'\n')
#cat('fac4',fac4,'\n')
ptot=p+d+s+dlam
#cat('ptot=',ptot,'\n')
#cat('prod.tmp',prod.tmp,'\n')
const=1
if (ptot > 0) {for(jp in 1:ptot) {const=const-phitot[jp]}}
#
# Calculate Box-Jenkins Forecasts
#
#
maconst=const*xbar
p1=max(ptot+1,q+1)
#
#Calculating residuals
#
#
for (i in p1:n) {resid[i]=x[i]
   if ( ptot > 0) {for (jp in 1:ptot) {resid[i]=resid[i]-phitot[jp]*x[i-jp]}}
   if (q > 0) {for (jq in 1:q) {resid[i]=resid[i]+theta[jq]*resid[i-jq]}}
                   resid[i]=resid[i]-maconst}
#
#Calculating forecasts
#
#
npn.ahead=n+n.ahead
xhat=rep(0,npn.ahead)
mm=n
if(lastn==TRUE) {mm=n-n.ahead}
for (i in 1:mm) {xhat[i]=x[i]}
for (h in 1:n.ahead) {
if (ptot > 0) {for (jp in 1:ptot) {xhat[mm+h]=xhat[mm+h]+phitot[jp]*xhat[mm+h-jp]}}
if (h<=q) {for(jq in h:q) {xhat[mm+h]=xhat[mm+h]-theta[jq]*resid[mm+h-jq]}}
                    xhat[mm+h]=xhat[mm+h]+maconst}
#
#
#   Calculate psi weights for forecasts limits
#
#
xi=psi.weights.wge(phitot,theta,lag.max=n.ahead)
#
#
#
#Setting up for plots
nap1=n.ahead+1
fplot=rep(0,nap1)
maxh=n+n.ahead
llplot=rep(0,nap1)
ulplot=rep(0,nap1)
f=rep(0,nap1)
ll=rep(0,nap1)
ul=rep(0,nap1)
wnv=0
xisq=rep(0,n.ahead)
se=rep(0,n.ahead)
se0=1
for (i in p1:n) {wnv=wnv+resid[i]**2}
wnv=wnv/(n-ptot)
xisq[1]=1
for (i in 2:n.ahead) {xisq[i]=xisq[i-1]+xi[i-1]^2}
for (i in 1:n.ahead) {se[i]=sqrt(wnv*xisq[i])}
fplot[1]=x[mm]
for (i in 1:n.ahead) {fplot[i+1]=xhat[mm+i]}
ulplot[1]=x[mm]
for (i in 1:n.ahead) { ulplot[i+1]=fplot[i+1]+1.96*se[i]}
llplot[1]=x[mm]
for (i in 1:n.ahead) { llplot[i+1]=fplot[i+1]-1.96*se[i]}
#
if(limits==FALSE) {
  if(lastn==TRUE) {maxp=max(x,xhat[1:n])
  minp=min(x,xhat[1:n])}
  else            {maxp=max(x,xhat)
  minp=min(x,xhat)}}
if(limits==TRUE) {minp=min(x,llplot)
maxp=max(x,ulplot)}

#cat(minp,maxp,'\n')
numrows <- 1
numcols <- 1
timelab <- 'Time'
valuelab <- ''
fig.width <- 5
fig.height <- 2.5
cex.labs <- c(.8,.7,.8)
#tiff(filename=filename,width=fig.width,height=fig.height,units='in',compression='none',res=350)
par(mfrow=c(numrows,numcols),mar=c(3.8,2.5,1,1))
t<-1:n;
np1=n+1
np.ahead=mm+n.ahead
tf<-mm:np.ahead
if (plot==TRUE) {plot(t,x,type='o',xaxt='n',yaxt='n',cex=.4,pch=16,cex.lab=.75,cex.axis=.75,lwd=.75,xlab='',ylab='',xlim=c(1,maxh),ylim=c(minp,maxp))
axis(side=1,cex.axis=.8,mgp=c(3,0.15,0),tcl=-.3);
axis(side=2,las=1,cex.axis=.8,mgp=c(3,.4,0),tcl=-.3)
mtext(side=c(1,2,1),cex=cex.labs,text=c(timelab,valuelab,""),line=c(.8,1.1,1.8))
points(tf,fplot,type='l',lty=3,cex=.4,lwd=2,pch=1);
if(limits==TRUE) {points(tf,ulplot,type='l',lty=3,cex=0.6,lwd=.75,pch=1)
points(tf,llplot,type='l',lty=3,cex=0.6,lwd=.75,pch=18) }
    }
np1=mm+1
nap1=n.ahead+1
f=fplot[2:nap1]
ll=llplot[2:nap1]
ul=ulplot[2:nap1]
out1=list(f=f,ll=ll,ul=ul,resid=resid,wnv=wnv,se=se,psi=xi,ptot=ptot,phitot=phitot)
return(out1)
}

