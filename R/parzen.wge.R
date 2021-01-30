
#
#  Program to calculate and plot the smoothed periodogram based on the Parzen window
#  
#
#
#
#
#
#
parzen.wge=function(x, dbcalc="TRUE", trunc=0, plot='TRUE')  
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
aut=acf(x,lag.max=n-1,plot=FALSE)$acf
#
#
# Periodogram
#
#
    freq=(1:floor(n/2))/n
    nf=length(freq)
    pgram=rep(0,nf)
    for (i in 1:nf)
  {
     cosvector=cos(2*pi*freq[i]*(1:(n-1)))
#cat('i, cosvector',i, cosvector,'\n')
       pgram[i]=aut[1]+2*sum(aut[-1]*cosvector)
#cat('i, freq[i],pgram[i]',i,freq[i],pgram[i],'\n')
  }
#
#
#list(freq=freq,pgram=pgram)
#} 
#
#
# Parzen Window
#
if (trunc == 0) {M=floor(2*sqrt(n))}
                else {M=trunc}
#
if(M>0) {
pzgram=rep(0,nf)
for (i in (1:nf))
    {cosvector=cos(2*pi*freq[i]*(1:(n-1)))
     weight=c( 
     c(1-6*((0:floor(M/2))/M)^2+6*((0:floor(M/2))/M)^3),
     c(2*(1-((floor(M/2)+1):M)/M)^3),
     rep(0,n-M-1) 
     )
     pzgram[i]=aut[1]*weight[1]+2*sum(aut[-1]*weight[-1]*cosvector)
     }
vlabel="dB"
if(dbcalc=="FALSE") {vlabel="Spectral Estimate (non-log)"}

if (dbcalc == "TRUE") {pzgram=10*log10(pzgram)}
if (plot=='TRUE') {
plot(freq,pzgram,type='l',xaxt='n',yaxt='n',pch=16,cex.lab=1,cex.axis=1,lwd=.75,xlab='',ylab='',
main=paste('Parzen Window Truncation point: M =',M))
axis(side=1,cex.axis=1,mgp=c(3,0.15,0),tcl=-.3);
axis(side=2,las=1,cex.axis=1,mgp=c(3,.4,0),tcl=-.3)
mtext(side=c(1,2,1),cex=cex.labs,text=c('Frequency',vlabel,''),line=c(1.5,1.3,2.1))
                 }
         }
#                      
out1=list(freq=freq,pzgram=pzgram)
return(out1)                       
               }                   
#