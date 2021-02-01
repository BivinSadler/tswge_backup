



#
artrans.wge=function(x,phi.tr,lag.max=25,plottr='TRUE') {
#
#   Perform an AR Transformation Y(t)=phi.tr(B)X(t) (see text, pages ...)
#   Examples:  To calculate Y(t)=(1-B)X(t)= X(t)-X(t-1) then phit.tr=c(1)
#              To calculate Y(t)=(1-1.6B+B^2)X(t) then phi.tr=c(1.6,-1)
#  
#
#  Input variables:
#  x=original realization
#  phi.tr=vector of transforming coefficients
#  lag.max=max lag for calulated sample autocorrelations
#  plottr=logical variable.  
#    If plottr=TRUE then ar.trans plots the original and transformed realizations and sample autocorrelations
#
# 
#  Output variables
#  y=transformed realization
#  aut.tr=sample autocorrelations of x for lags 0 to lag.max
#
# Calculate sample autocorrelations for x
#
p.tr=length(phi.tr)
n=length(x)
nmp=n-p.tr
nmp1=nmp-1
if(lag.max > nmp1) {lag.max=nmp1}
aut=acf(x,lag.max,plot=FALSE)$acf
lag.max1=lag.max+1
k=0:lag.max
aut.orig=aut[1:lag.max1]
#
# Calculate transformed realization y
#
y=rep(0,nmp)
#cat(p.tr,nmp,y,'\n')
for (t in p.tr:n) {
y[t-p.tr]=x[t]
for  (i in 1:p.tr)       {
y[t-p.tr]=y[t-p.tr]-phi.tr[i]*x[t-i]
                   }
                  }
#
#  Calculate sample autocorrelations for y
#
aut=acf(y,lag.max,plot=FALSE)$acf
aut.tr=aut[1:lag.max1]
#
  if (plottr == TRUE) {

# Preparing to plot the 4 plots
#
#
#
numrows <- 2
numcols <- 2
fig.width <- 5.5
fig.height <- 4.5
cex.labs <- c(.8,.7,.9)
par(mfrow=c(numrows,numcols),mar=c(3.8,2.5,1,1))
#
#
#  Plot X
#
tt=1:n
plot(tt,x,type='o',xaxt='n',yaxt='n',cex=0.5,pch=16,cex.lab=.75,cex.axis=.75,lwd=.75,xlab='',ylab='')
axis(side=1,cex.axis=.9,mgp=c(3,0.15,0),tcl=-.3);
axis(side=2,las=1,cex.axis=.9,mgp=c(3,.4,0),tcl=-.3)
mtext(side=c(1,2,1),cex=cex.labs,text=c('Time','','Realization'),line=c(1,1.1,2.1))
#                 }
#
#
# Plot Sample Autocorrelations for X
#
#
#aut.orig=aut[1:lag.max1]
plot(k,aut.orig,type='h',xaxt='n',yaxt='n',cex=0.0,cex.lab=.75,cex.axis=.75,lwd=.75,xlab='',ylab='',ylim=c(-1,1))
abline(h=0)
axis(side=1,cex.axis=.8,mgp=c(3,0.15,0),tcl=-.3);
axis(side=2,las=1,cex.axis=.8,mgp=c(3,.4,0),tcl=-.3)
mtext(side=c(1,2,1),cex=cex.labs,text=c('Lag','','Sample Autocorrelations'),line=c(1,1.1,2.1))    
#
#
#  Plot y=phi.tr(X)
if (plottr=='TRUE') {
#
tt=1:nmp
plot(tt,y,type='o',xaxt='n',yaxt='n',cex=0.5,pch=16,cex.lab=.75,cex.axis=.75,lwd=.75,xlab='',ylab='')
axis(side=1,cex.axis=.9,mgp=c(3,0.15,0),tcl=-.3);
axis(side=2,las=1,cex.axis=.9,mgp=c(3,.4,0),tcl=-.3)
mtext(side=c(1,2,1),cex=cex.labs,text=c('Time','','Transformed Realization'),line=c(1,1.1,2.1))
#                 }
#
#
# Plot Sample Autocorrelations for Y
#
#
k=0:lag.max
plot(k,aut.tr,type='h',xaxt='n',yaxt='n',cex=0.0,cex.lab=.75,cex.axis=.75,lwd=.75,xlab='',ylab='',ylim=c(-1,1))
abline(h=0)
axis(side=1,cex.axis=.8,mgp=c(3,0.15,0),tcl=-.3);
axis(side=2,las=1,cex.axis=.8,mgp=c(3,.4,0),tcl=-.3)
mtext(side=c(1,2,1),cex=cex.labs,text=c('Lag','','Sample Autocorrelations'),line=c(1,1.1,2.1))   
                              }
}
return(y) 
}
#
#
#
#
