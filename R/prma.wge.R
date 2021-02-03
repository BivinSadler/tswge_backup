prma.wge=function(x,order=3,k.ahead=1)
{
n=length(x)
nk=n+k.ahead
n1=n+1
ord1=order+1
t1=ord1:nk
t=1:n
ma=rep(0,nk)
sum=0
for(i in 1:order) sum=sum+x[order-i+1]
ma[ord1]=sum/order
ord2=order+2
for (i in ord2:n1) ma[i]=ma[i-1]-x[i-order-1]/order+x[i-1]/order
if(k.ahead>1) {
x[n1]=ma[n1]
for(i in 2:k.ahead) {ma[n1+i-1]=ma[n1+i-2]-x[n1+i-1-order-1]/order+x[n1+i-2]/order
x[n1+i-1]=ma[n1+i-1]
}
}
fig.width <- 7.5
fig.height <-3.8
par(mfrow=c(2,2),mar=c(3,2.2,1,.5))
cex.labs <- c(1.5,1.5,1.5)
plot(t,x[1:n],type='o',xaxt='n',yaxt='n',cex=0.5,pch=16,cex.lab=.75,lwd=.75,xlab='',ylab='',xlim=c(0,nk))
points(t1,ma[ord1:nk],type='o',cex=.8,pch=15,lwd=2)
axis(side=1,cex.axis=1.4,mgp=c(3,0.45,0),tcl=-.3)
axis(side=2,las=1,cex.axis=1.4,mgp=c(3,.45,0),tcl=-.3)
mtext(side=c(1,2,1),cex=cex.labs,text=c('Time','','Data and MA Predictors'),line=c(1.7,1.8,3.5))
#
out1=list(x=x[1:n],ma=ma,k.ahead=k.ahead)
return(out1)       
}
##
