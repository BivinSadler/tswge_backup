#  Program to calculate and plot the periodogram
#  
#
#
#
#
#
#
period.wge=function(x,dbcalc='TRUE',plot='TRUE')
{
#
# x=vector containing realization to be analyzed
# dbcalc is a logical variable, dbcalc="TRUE" (default) #specifies that the periodogram will be calculated on the log #scale (in dB), 
#    i.e. 10log10(periodogram)
# plot is a logical variable, plot="TRUE" (default) specifies #that the periodogram will be plotted 
#    (either log or non-log as specified by db)
#

#fig.width <- 5.5
#fig.height <- 4.5
cex.labs <- c(1,1,1)
#par(mfrow=c(numrows,numcols),mar=c(3.8,2.5,2.5,1))
n=length(x)
t=1:n
aut=acf(x,lag.max=n-1,plot=FALSE)$acf
#
#
# Calculate Periodogram
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

vlabel='dB' 
if (dbcalc=='FALSE') {vlabel='Periodogram (non-log)'}
if (dbcalc == 'TRUE') {pgram=10*log10(pgram)}
if(plot=="TRUE") {
nf=length(pgram)
min.per=min(pgram[1:nf])
plot(freq,pgram,type='n',xaxt='n',yaxt='n',cex=0.5,pch=16,cex.lab=.75,cex.axis=.75,lwd=.75,xlab='',ylab='', main='Periodogram')
axis(side=1,cex.axis=1,mgp=c(3,0.15,0),tcl=-.3);
axis(side=2,las=1,cex.axis=1,mgp=c(3,.4,0),tcl=-.3)
mtext(side=c(1,2,1),cex=cex.labs,text=c('Frequency',vlabel,''),line=c(1.5,1.3,2.1))
for (i in 1:nf) {segments(freq[i],min.per,freq[i],pgram[i])
                }   
              }     
#
out1=list(freq=freq,pgram=pgram)
return(out1) 
               }                   

#