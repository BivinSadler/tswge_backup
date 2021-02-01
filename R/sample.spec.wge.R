
#
#  Program to calculate and plot the smoothed periodogram based on the Parzen window
#  
#
#
#
#
#
#
sample.spec.wge=function(x, dbcalc="TRUE", plot='TRUE')  
{
#
# x=vector containing realization to be analyzed
# dbcalc is a logical variable, dbcalc="TRUE" (default) specifies that the smoothed periodogram will be calculated on the log scale 
#    (in dB), i.e. 10log10(periodogram)
# trunc specifies the truncation point. if M=0 (default) the the calculation is based on the default truncation point M=2(n)^.5
# trunc>0 is a user-specified truncation point
# plot="TRUE" (default) specifies that the smoothed periodogram will be plotted 
#    (either log or non-log as specified by dbcalc)
#
#numcols <- 2
#fig.width <- 5.5
#fig.height <- 4.5
cex.labs <- c(1,1,1)
#par(mfrow=c(numrows,numcols),mar=c(3.8,2.5,2.5,1))
n=length(x)
t=1:n
nm1=n-1
aut=acf(x,lag.max=nm1,plot=FALSE)$acf
#
#
# Sample Spectral Density
#
#
f=1:251
freq=(f-1)/500
nf=251
sample.sp=rep(1,251)
    for (i in 2:nf)
    {
       for(k in 1:nm1) {sample.sp[i]=sample.sp[i]+2*aut[k+1]*cos(2*pi*freq[i]*k)}
}
vlabel="dB"
if(dbcalc=="FALSE") {vlabel="Spectral Estimate (non-log)"}

if (dbcalc == "TRUE") {sample.sp=10*log10(sample.sp)}
if (plot=='TRUE') {
plot(freq,sample.sp,type='l',xaxt='n',yaxt='n',pch=16,cex.lab=1,cex.axis=1,lwd=.75,xlab='',ylab='')
axis(side=1,cex.axis=1,mgp=c(3,0.15,0),tcl=-.3);
axis(side=2,las=1,cex.axis=1,mgp=c(3,.4,0),tcl=-.3)
mtext(side=c(1,2,1),cex=cex.labs,text=
c('Frequency',vlabel,''),line=c(1.5,1.3,2.1))
                 }
#                      
out1=list(freq=freq,sample.sp=sample.sp)
return(out1)                       
               }                   


#