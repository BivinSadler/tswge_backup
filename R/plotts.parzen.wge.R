
#
#
#  Program to calculate and plot the periodogram and up to three smoothed spectral estimators (either using the Parzen window)
#  
#
#
#
#
#
#
plotts.parzen.wge=function(x,m2=c(0,0))
{
#
# x=vector containing realization to be analyzed
# window.type can take on the values "Parzen" (default) or "Bartlett" 
# m2 is a 2-component vector of alternative, user supplied, truncation points.
#    Example: m2=c(10,24) indicates that in addition to the default truncation point, the smoothed spectral estimator
#                         is to be calculated using truncation points 10 and 24
#             m2=c(0,0) indicates that no additional truncation points are to be used
#             m2=c(10,0) indicates the use of one additional truncation point (10) 
#
#
#  Output
#
#  For the command spec=spec.sample.atsa(x,M2)
# 
#  spec$freq - vector containing the frequencies at which the spectrum is calculated (1/n, 2/n, ... .5 (n even)) (1/n, 2/n, ... (n-1)/2n (n odd))
#               (the frequencies in spec$freq are the harmonics of the fundamental frequency 1/n)
#  spec$db  - vector containing the periodogram (in db) calulated at the frequencies in spec$freq
#  spec$dbz - vector containing the Parzen spectral density (in db) calculated at the frequencies in spec$freq with m=2(n).5
#  spec$dbz1 - (optionally) vector containing the Parzen spectral density (in db) calculated at the frequencies in spec$freq at m2[1]
#  spec$dbz2 - (optionally) vector containing the Parzen spectral density (in db) calculated at the frequencies in spec$freq at m2[2]
numrows <- 2
numcols <- 2
#fig.width <- 5.5
#fig.height <- 5.5
cex.labs <- c(.8,.8,.8)
par(mfrow=c(numrows,numcols),mar=c(2.8,2.5,2.5,1))
n=length(x)
t=1:n
aut=acf(x,lag.max=n-1,plot=FALSE)$acf
#
#
# Plot Periodogram
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
 
db=10*log10(pgram)
nf=length(db)
min.per=min(db[1:nf])
plot(freq,db,type='n',xaxt='n',yaxt='n',cex=0.8,pch=16,cex.lab=.75,cex.axis=.75,lwd=.75,xlab='',ylab='', main='Periodogram',cex.main=1.2)
axis(side=1,cex.axis=1,mgp=c(3,0.15,0),tcl=-.3);
axis(side=2,las=1,cex.axis=1,mgp=c(3,.4,0),tcl=-.3)
mtext(side=c(1,2,1),cex=cex.labs,text=c('Frequency','dB',''),line=c(1.5,1.3,2.1))
for (i in 1:nf) {segments(freq[i],min.per,freq[i],db[i])
                 }   
#
#
# Parzen Window
#

#
#
for (ii in 1:3) {
     if(ii == 1) {m=floor(2*sqrt(n))}
                 else {ii1=ii-1
                       m=m2[ii1]}

if(m>0) {
pzgram=rep(0,nf)
for (i in (1:nf))
    {cosvector=cos(2*pi*freq[i]*(1:(n-1)))
     weight=c( 
     c(1-6*((0:floor(m/2))/m)^2+6*((0:floor(m/2))/m)^3),
     c(2*(1-((floor(m/2)+1):m)/m)^3),
     rep(0,n-m-1) 
     )
     pzgram[i]=aut[1]*weight[1]+2*sum(aut[-1]*weight[-1]*cosvector)
     }
dbz=10*log10(pzgram)
plot(freq,dbz,type='l',xaxt='n',yaxt='n',cex=0.8,pch=16,cex.lab=.75,cex.axis=.75,lwd=.75,xlab='',ylab='',main='')
title(main=paste('Parzen Window Truncation point: M =',m),cex.main=1.2)
axis(side=1,cex.axis=1,mgp=c(3,0.15,0),tcl=-.3);
axis(side=2,las=1,cex.axis=1,mgp=c(3,.4,0),tcl=-.3)
mtext(side=c(1,2,1),cex=cex.labs,text=c('Frequency','dB',''),line=c(1.5,1.3,2.1))

         }
#                      
if(ii==1){dbz=dbz}
if(ii==2){dbz1=dbz}
if(ii==3){dbz2=dbz}
                 }
out1=list(freq=freq,db=db,dbz=dbz,dbz1=dbz1,dbz2=dbz2)
return(out1)                       
               }                   
#