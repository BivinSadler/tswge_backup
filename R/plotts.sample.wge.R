
###be cautious whether I have loaded tSA before calling acf - they give different results (regarding counting lag #zero, etc.)###
#
# 
#
#
# function plotts.sample.wge
#
# Plot a time series realization, x, sample autocorrelations, periodogram, and Parzen window smoother (based on #default truncation point 2(n)^.5)
#
#
plotts.sample.wge=function(x,lag.max=25,trunc=0,arlimits=FALSE,speclimits=c(0,0),periodogram=FALSE) 
{
#
#
#     x is a vector of length n containing the time series realization  
#     lag.max  is the maximum lag at which to calculate the sample autocorrelations 
#     
#     trunc (integer =0) specifies the truncation point.
#          trunc=0 (default) (or no value for trunc is specified) indicates that the 
#                   default truncation point is  M=2(n)^.5 
#          trunc>0 is a user specified truncation point, i.e. calculations will be based on 
#          ar limits is a logical variable specifying whether 95% limit lines will be included on sample autocorrelation plots
#
xbar=mean(x)
cex.labs=c(1,1,1)
if(periodogram == TRUE) {
dev.new(width=8,height=6.5)
numrows=2
numcols=2
par(mfrow=c(numrows,numcols),mar=c(4.3,3.5,1.6,1))}
#
if(periodogram == FALSE) {
dev.new(width=11,height=3)
numrows=1
numcols=3
par(mfrow=c(numrows,numcols),mar=c(4.3,3.5,1.6,1))}
#
#
n=length(x)
nm1=n-1
if(lag.max > nm1) {lag.max=nm1}
t=1:n
naut=lag.max*1
naut1=naut+1
k=0:naut
aut=acf(x,lag.max=n-1,plot=FALSE)$acf
#aut=acf(x,lag.max=naut,plot=FALSE) 
#
#
#PLOT DATA
#
#
if (n < 200) {plot(t,x,type='o',xaxt='n',yaxt='n',cex=0.5,pch=16,cex.lab=.75,lwd=.75,xlab='',ylab='')}
if (n >= 200) {plot(t,x,type='l',xaxt='n',yaxt='n',cex=0.5,pch=16,cex.lab=.75,lwd=.75,xlab='',ylab='')}
if(periodogram==FALSE) {axis(side=1,cex.axis=1.4,mgp=c(3,0.45,0),tcl=-.3)};
if(periodogram==TRUE) {axis(side=1,cex.axis=1.2,mgp=c(3,0.45,0),tcl=-.3)};
if(periodogram==FALSE) {axis(side=2,las=1,cex.axis=1.4,mgp=c(3,.45,0),tcl=-.3)}
if(periodogram==TRUE) {axis(side=2,las=1,cex.axis=1.2,mgp=c(3,.45,0),tcl=-.3)}
mtext(side=c(1,2,1),cex=cex.labs,text=c('Time','','Realization'),line=c(1.7,1.8,3))
#                 }
#
#
# SAMPLE AUTOCORRELATIOS
#
#
ul=2/sqrt(n)
ll=-ul
autplt=aut[1:naut1]
plot(k,autplt,type='h',xaxt='n',yaxt='n',cex=0.0,cex.lab=.75,cex.axis=.75,lwd=.75,xlab='',ylab='',ylim=c(-1,1))
abline(h=0)
if(arlimits==TRUE) {abline(h=ul,lty=2)
                    abline(h=-ul,lty=2)}
if(periodogram==FALSE) {axis(side=1,cex.axis=1.4,mgp=c(3,0.45,0),tcl=-.3)};
if(periodogram==TRUE) {axis(side=1,cex.axis=1.2,mgp=c(3,0.45,0),tcl=-.3)};
if(periodogram==FALSE) {axis(side=2,las=1,cex.axis=1.4,mgp=c(3,.45,0),tcl=-.3)}
if(periodogram==TRUE) {axis(side=2,las=1,cex.axis=1.2,mgp=c(3,.45,0),tcl=-.3)}
mtext(side=c(1,2,1),cex=cex.labs,text=c('Lag','','Sample Autocorrelations'),line=c(1.7,1.8,3))
#                }
#

#
#
# PARZEN WINDOW
#
#
  freq=(1:floor(n/2))/n
  nf=length(freq)
#
if (trunc == 0) {M=floor(2*sqrt(n))}
                else {M=trunc}
#
#
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
dbz=10*log10(pzgram)
if(sum(speclimits^2)==0) speclimits=c(min(dbz),max(dbz))
plot(freq,dbz,type='l',xaxt='n',yaxt='n',cex=0.5,pch=16,cex.lab=.75,cex.axis=.75,lwd=.75,xlab='',ylab='',ylim=speclimits)
if(periodogram==FALSE) {axis(side=1,cex.axis=1.4,mgp=c(3,0.45,0),tcl=-.3)};
if(periodogram==TRUE) {axis(side=1,cex.axis=1.2,mgp=c(3,0.45,0),tcl=-.3)};
if(periodogram==FALSE) {axis(side=2,las=1,cex.axis=1.4,mgp=c(3,.45,0),tcl=-.3)}
if(periodogram==TRUE) {axis(side=2,las=1,cex.axis=1.2,mgp=c(3,.45,0),tcl=-.3)}
mtext(side=c(1,2,1),cex=cex.labs,text=c('Frequency','dB','Parzen Window'),line=c(1.7,1.8,3))
#               }
#
#
if(periodogram == TRUE) {
# PLOT PERIODOGRAM
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
db=10*log10(pgram)
nf=length(db)
min.per=min(db[1:nf])
if(n <= 200) {plot(freq,db,type='n',xaxt='n',yaxt='n',cex=0.5,pch=16,cex.lab=.75,cex.axis=.75,lwd=.75,xlab='',ylab='',ylim=speclimits)}
if(n > 200) {plot(freq,db,type='l',xaxt='n',yaxt='n',cex=0.5,pch=16,cex.lab=.75,cex.axis=.75,lwd=.75,xlab='',ylab='',ylim=speclimits)}
if(periodogram==FALSE) {axis(side=1,cex.axis=1.4,mgp=c(3,0.45,0),tcl=-.3)};
if(periodogram==TRUE) {axis(side=1,cex.axis=1.2,mgp=c(3,0.45,0),tcl=-.3)};
if(periodogram==FALSE) {axis(side=2,las=1,cex.axis=1.4,mgp=c(3,.45,0),tcl=-.3)}
if(periodogram==TRUE) {axis(side=2,las=1,cex.axis=1.2,mgp=c(3,.45,0),tcl=-.3)}
mtext(side=c(1,2,1),cex=cex.labs,text=c('Frequency','dB','Periodogram'),line=c(1.7,1.8,3))
for (i in 1:nf) {segments(freq[i],min.per,freq[i],db[i])
                }   
}
#
#
#
#
if (periodogram == FALSE) {out1=list(xbar=xbar,autplt=autplt,freq=freq,dbz=dbz)}
if (periodogram == TRUE) {out1=list(xbar=xbar,autplt=autplt,freq=freq,dbz=dbz,db=db)}
return(out1)                       
}

#
#

