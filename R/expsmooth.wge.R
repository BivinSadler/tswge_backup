expsmooth.wge=function(x,alpha=NULL,n.ahead=0,plot=TRUE)
#registerS3method("exp","default",base::exp)
#registerS3method("exp","smooth.wge",expsmooth.wge)
{
alpha=alpha
n=length(x)
npn.ahead=n+n.ahead
xhat=rep(0,npn.ahead)
#
# Calculating 1-step ahead Forecasts
#
#
npn.ahead=n+n.ahead
xhat=rep(0,npn.ahead)
x.hw=HoltWinters(x,alpha=alpha,beta=FALSE,gamma=FALSE)
u=rep(0,length(x))
for(j in 2:length(x)) u[j]=x.hw$fitted[,1][j-1]
u[1]=x[1]
numrows <- 1
numcols <- 1
fig.width <- 5
fig.height <- 1.5
cex.labs <- c(1.2,1.2,1.2)
par(mfrow=c(numrows,numcols),mar=c(3.8,3.5,1,1))
#
if(n.ahead==0) {
t=1:n
plot(t,x,type='o',xaxt='n',yaxt='n',cex=.8,pch=16,cex.lab=1,cex.axis=1,lwd=1,xlab='',ylab='',col=1)
points(t,u,type='l',pch=2,lwd=2,cex=.8)
axis(side=1,cex.axis=1.1,mgp=c(3,0.15,0),tcl=-.3);
axis(side=2,las=1,cex.axis=1.1,mgp=c(3,.4,0),tcl=-.3)
}
#
if(n.ahead>0){
f=predict(x.hw,n.ahead=n.ahead)
for(j in 1:n.ahead) u[n+j]=f[j]
plot(x,type='o',xaxt='n',yaxt='n',cex=.8,pch=16,cex.lab=1,cex.axis=1,lwd=1,xlab='',ylab='',xlim=c(1,npn.ahead),col=1)
points(u,type='o',pch=2,lwd=.8,cex=.8)
axis(side=1,cex.axis=1.1,mgp=c(3,0.15,0),tcl=-.3);
axis(side=2,las=1,cex.axis=1.1,mgp=c(3,.4,0),tcl=-.3)
}
out1=list(alpha=x.hw$alpha,u=u)
}
