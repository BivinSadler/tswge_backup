fore.aruma.wge=function(x,phi=0,theta=0,d=0,s=0,lambda=0,n.ahead=5,lastn=FALSE,plot=TRUE,alpha=.05,limits=TRUE)
{
#x=ngen.aruma.wge(n=100,phi=c(0,-.9),d=1,s=4,sn=50)
#n=100
#phi=c(.0534,-.9291)
#theta=0
#lambda=0
#d=1
#s=4
#plot=TRUE
#limits=FALSE
#n.ahead=20
n=length(x)
p=length(phi)
if(p == 1)
   {if(phi==0) {p=0
            fac1=0}}
if(sum(phi^2)>0) {fac1=phi}
q=length(theta)
if(sum(theta^2)==0) {q=0}
npn.ahead=n+n.ahead
xhat=rep(0,npn.ahead)
xbar=mean(x)
if (d > 3) {cat('d>3','\n')
            return(d)}
if (d==0) {diffac=0}
if (d==1) {diffac=1}
if (d==2) {diffac=c(2,-1)}
if (d==3) {diffac=c(3,-3,1)}
if (d<4) {fac2=diffac} 
cat('s=',s,'\n') 
if(s==0) {fac3=0}
if(s > 0) {seas=rep(0,s)
           seas[s]=1
           fac3=seas}
#cat('seas',seas,'\n')
dlam=length(lambda)
#cat('lambda',lambda,'\n')
if(sum(lambda^2)==0) {dlam=0
            fac4=0}
if(sum(lambda^2)>0) {fac4=lambda}
#cat('lambda',lambda,'\n')
prod.fore.tmp=mult.wge(fac1,fac2,fac3,fac4)
prod.res.tmp=mult.wge(fac2,fac3,fac4)
phitot.fore=prod.fore.tmp$model.coef
phitot.res=prod.res.tmp$model.coef
#cat('fac1',fac1,'\n')
#cat('fac2',fac2,'\n')
#cat('fac3',fac3,'\n')
#cat('fac4',fac4,'\n')
ptot.fore=p+d+s+dlam
ptot.res=d+s+dlam
cat('ptot.res=',ptot.res,'\n')
cat('phitot.res',phitot.res,'\n')
cat('ptot.fore=',ptot.fore,'\n')
cat('phitot.fore',phitot.fore,'\n')
y.arima=x
lns=length(phitot.res)
if(lns>0) {y.arma=artrans.wge(x,phi.tr=phitot.res)}
#
#BACKCAST RESIDUALS BASED ON FIT TO TRANSFORMED ARMA DATA
#

res=rep(0,n)
res
res=backcast.wge(y.arma,phi=phi,theta=theta,n.back=50)
res
ntot.res=length(y.arma)
resid=rep(0,n)
resid
for (i in 1:ntot.res){resid[ptot.res+i]=res[i]}
resid
#
#
# Calculate Box-Jenkins Forecasts
const=1
if (ptot.fore> 0) {for(jp in 1:ptot.fore) {const=const-phitot.fore[jp]}}
#
#
maconst=const*xbar
p1=max(ptot.fore+1,q+1)
#
#Caclulating residuals (old)
#
#
#for (i in p1:n) {resid[i]=x[i]
#   if ( ptot > 0) {for (jp in 1:ptot) {resid[i]=resid[i]-phitot[jp]*x[i-jp]}}
#   if (q > 0) {for (jq in 1:q) {resid[i]=resid[i]+theta[jq]*resid[i-jq]}}
#                   resid[i]=resid[i]-maconst}
#
#CALCULATING FORECASTS
#
#
npn.ahead=n+n.ahead
xhat=rep(0,npn.ahead)
mm=n
if(lastn==TRUE) {mm=n-n.ahead}
for (i in 1:mm) {xhat[i]=x[i]}
for (h in 1:n.ahead) {
if (ptot.fore > 0) {for (jp in 1:ptot.fore) {xhat[mm+h]=xhat[mm+h]+phitot.fore[jp]*xhat[mm+h-jp]}}
if (h<=q) {for(jq in h:q) {xhat[mm+h]=xhat[mm+h]-theta[jq]*resid[mm+h-jq]}}
                    xhat[mm+h]=xhat[mm+h]+mean(x[1:mm])*const}
#
#
#   Calculate psi weights for forecasts limits
#
#
xi=psi.weights.wge(phitot.fore,theta,lag.max=n.ahead)
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
resid
for (i in p1:n) {wnv=wnv+resid[i]**2}
wnv=wnv/(n-ptot.res)
xisq[1]=1
for (i in 2:n.ahead) {xisq[i]=xisq[i-1]+xi[i-1]^2}
for (i in 1:n.ahead) {se[i]=sqrt(wnv*xisq[i])}
fplot[1]=x[mm]
for (i in 1:n.ahead) {fplot[i+1]=xhat[mm+i]}
ulplot[1]=x[mm]
for (i in 1:n.ahead) { ulplot[i+1]=fplot[i+1]-qnorm(alpha/2)*se[i]}
llplot[1]=x[mm]
for (i in 1:n.ahead) { llplot[i+1]=fplot[i+1]+qnorm(alpha/2)*se[i]}
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
cex.labs <- c(1.2,1.2,1.2)
par(mfrow=c(numrows,numcols),mar=c(9,4,3,2))
t<-1:n;
np1=n+1
np.ahead=mm+n.ahead
tf<-mm:np.ahead
if (plot==TRUE) {plot(t,x,type='o',xaxt='n',yaxt='n',cex=.8,pch=16,cex.lab=1,cex.axis=1,lwd=.75,xlab='',ylab='',xlim=c(1,maxh),ylim=c(minp,maxp))
axis(side=1,cex.axis=1.1,mgp=c(3,0.15,0),tcl=-.3);
axis(side=2,las=1,cex.axis=1.1,mgp=c(3,.4,0),tcl=-.3)
mtext(side=c(1,2,1),cex=cex.labs,text=c(timelab,valuelab,""),line=c(1.2,2.1,1.8))
points(tf,fplot,type='o',lty=1,cex=.8,lwd=.5,pch=1);
if(limits==TRUE) {points(tf,ulplot,type='l',lty=3,cex=0.6,lwd=.75,pch=1)
points(tf,llplot,type='l',lty=3,cex=0.6,lwd=.75,pch=18) }
    }
np1=mm+1
nap1=n.ahead+1
f=fplot[2:nap1]
ll=llplot[2:nap1]
ul=ulplot[2:nap1]
out1=list(f=f,ll=ll,ul=ul,resid=resid,wnv=wnv,xbar=xbar,se=se,psi=xi,ptot.fore=ptot.fore,phitot.fore=phitot.fore)
return(out1)
}

