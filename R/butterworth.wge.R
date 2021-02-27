#
# Butterworth Filter
#
# Use package Signal using filtfilt after calling butter
#
#
butterworth.wge=function(x,order=4,type,cutoff,plot=TRUE) {
#  x is time series realization
#  order is order of Butterworth filter
#  type= either "high", "low", "stop", or "band"
#  cutoff
#   if type="high" or "low" the cutoff is the scalar cutoff frequency
#   if type="stop" or "band" then cutoff is a 2-component vector containing thelower and upper cutoff frequencies
#plot
#  if plot=TRUE then plots of original and reinterpolated series are plotted
#
#
n=length(x)
cutoff2=cutoff*2
bb=signal::butter(order,cutoff2,type)
x.filt=signal::filtfilt(bb$b,bb$a,x)
#                                    
if(plot==TRUE) {
max.x=max(x[1:n])
min.x=min(x[1:n])
max.f=max(x.filt[1:n])
min.f=min(x.filt[1:n])
max.both=max(max.x,max.f)
min.both=min(min.x,min.f)
par(mfrow=c(2,1),mar=c(3,2.5,1,.5))
cex.labs <- c(.8,.7,.8)
t=1:n
fig.width=6
fig.height=3
#
plot(t,x,type='l',xaxt='n',yaxt='n',cex=0.4,cex.lab=.75,cex.axis=.75,lwd=.75,xlab='',ylab='',ylim=c(min.both,max.both))
axis(side=1,cex.axis=.8,mgp=c(3,0.15,0),tcl=-.3);
axis(side=2,las=1,cex.axis=.8,mgp=c(3,.4,0),tcl=-.3)
mtext(side=c(1,2,1),cex=cex.labs,text=c('Time','', 'Original series'),line=c(1,1.1,2))
#
plot(t,x.filt,type='l',xaxt='n',yaxt='n',cex=0.4,cex.lab=.75,cex.axis=.75,lwd=.75,xlab='',ylab='',ylim=c(min.both,max.both))
axis(side=1,cex.axis=.8,mgp=c(3,0.15,0),tcl=-.3);
axis(side=2,las=1,cex.axis=.8,mgp=c(3,.4,0),tcl=-.3)
mtext(side=c(1,2,1),cex=cex.labs,text=c('Time','','Filtered Data'),line=c(1,1.1,2))
}
#
#
out1=list(x.filt=x.filt)
return(out1)                                             }
#
#
#
