true.farma.aut.wge=function(d=0,phi=0,theta=0,lag.max=50,trunc=1000,vara=1,plot=TRUE)
{
ar=phi 
ma=theta
maxlag=lag.max
acv_temp=true.arma.aut.wge(ar,ma,trunc,vara,plot=FALSE)
acv_V=acv_temp$acv
acv_X=c()
 acv_X[1]=gamma(1-2*d)/(gamma(1-d))^2          
 for (i in 2:(maxlag+trunc+1))
 { 
  acv_X[i]=acv_X[i-1]*(i-1-1+d)/(i-1-d)
 }

 acv_Y=c()
 for (k in 1:(maxlag+1))
 {
  xx=c()
  for (j in 1:(trunc+1))
  { 
  xx[j]=acv_V[j]*(acv_X[j+k-1]+acv_X[abs(k-j)+1])
  }

  acv_Y[k]=sum(xx)-acv_V[1]*acv_X[k]
 }
acf=acv_Y/acv_Y[1]
#
# plot true autocorrelations
#
if(plot=="TRUE") {k=0:maxlag
numrows <- 4
numcols <- 1
fig.width <- 4
fig.height <- 3
par(mfrow=c(numrows,numcols),mar=c(2.5,2.5,1,1.5))
cex.labs <- c(.9,.8,.9)
plot(k,acf,type='h',xaxt='n',yaxt='n',cex=0.4,cex.lab=.85,cex.axis=.85,lwd=.85,xlab='',ylab='',ylim=c(-1,1))
abline(h=0)
axis(side=1,cex.axis=.85,mgp=c(3,0.15,0),tcl=-.3);
axis(side=2,las=1,cex.axis=.85,mgp=c(3,.4,0),tcl=-.3)
mtext(side=c(1,2,1),cex=cex.labs,text=c('Lag','','True Autocorrelations'),line=c(1,1.1,2.1))
}
out1=list(acf=acf,acv=acv_Y)
return(out1)
}



