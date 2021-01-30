library(MASS)
true.garma.aut.wge=function(u,lambda,phi=0,theta=0,lag.max=50,vara=1,plot=TRUE)
{
 if(sum(abs(phi))!=0) {p=length(phi)} else p=0
 if(sum(abs(theta))!=0) {q=length(theta)} else q=0

 acv=c()
 for (i in 1:(lag.max+1))
 { 
  newfunc=function(f)
  {
   if(q==0) {comp1=1} else comp1=Mod(sum(c(1,-theta)*(exp(2*pi*(1i)*f)^c(0:q))))^2
   if(p==0) {comp2=1} else comp2=Mod(sum(c(1,-phi)*(exp(2*pi*(1i)*f)^c(0:p))))^2 
   comp3=prod(  ( 4*(cos(2*pi*f)-u)^2 )^(-lambda)  )
   if(f==acos(u)/(2*pi)) {result1=0}
   else result1=vara*comp1/comp2*comp3
   result.temp=result1*cos(2*pi*f*(i-1))
   return(result.temp)
  }
  acv[i]=2*MASS::area(newfunc,0,.5) 
 }  
acf=acv/acv[1] 
#
#
# plot true autocorrelations
#
if(plot=="TRUE") {k=0:lag.max
numrows <- 4
numcols <- 1
fig.width <- 4
fig.height =3
par(mfrow=c(numrows,numcols),mar=c(2.5,2.5,1,1.5))
cex.labs <- c(.9,.8,.9)
plot(k,acf,type='h',xaxt='n',yaxt='n',cex=0.4,cex.lab=.75,cex.axis=.75,lwd=.75,xlab='',ylab='',ylim=c(-1,1))
abline(h=0)
axis(side=1,cex.axis=.8,mgp=c(3,0.15,0),tcl=-.3);
axis(side=2,las=1,cex.axis=.8,mgp=c(3,.4,0),tcl=-.3)
mtext(side=c(1,2,1),cex=cex.labs,text=c('Lag','','True Autocorrelations'),line=c(1,1.1,2.1))
}
#
out1=list(acf=acf,acv=acv)
return(out1)
}
