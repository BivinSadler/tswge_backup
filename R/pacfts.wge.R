
###be cautious whether I have loaded tSA before calling acf - they give different results (regarding counting lag #zero, etc.)###
#
# 
#
#
# function pacf.wge
#
# Plot a pacf of time series realization
#
#
pacfts.wge=function(x,lag.max=5,plot=TRUE,na.action,limits=FALSE, method='yw')
{
#
#
#     x is a vector of length n containing the time series realization  
#     lag.max  is the maximum lag at which to calculate the partial autocorrelations 
#     limits is a logical variable specifying whether 95% limit lines will be included on sample autocorrelation plots
#

cex.labs <- c(1.4,1.4,1.4)
#
numrows <- 1
numcols <- 1
par(mfrow=c(numrows,numcols),mar=c(4.2,4.2,8,1))
#
n=length(x)
nm1=n-1
pacf_burg=rep(0,nm1)
pacf_yw=rep(0,nm1)
pacf_mle=rep(0,nm1)
if(lag.max > nm1) {lag.max=nm1}
t=1:n
naut=lag.max
k=1:naut
if (method=="yw") {
   for(order in 1:lag.max) {
      temp=est.ar.wge(x,p=order,method='yw',factor=FALSE)
      pacf_yw[order]=temp$phi[order]
                            }
                   }
if (method=="burg") {
   for(order in 1:lag.max) {
      temp=est.ar.wge(x,p=order,method='burg',factor=FALSE)
      pacf_burg[order]=temp$phi[order]
                            }
                   }
if (method=="mle") {
   for(order in 1:lag.max) {
      temp=est.ar.wge(x,p=order,method='mle',factor=FALSE)
      pacf_mle[order]=temp$phi[order]
                            }
                   }

#
#
# Plot Partial Autocorrelations
#
#
ul=2/sqrt(n)
ll=-ul
#
if(method=='yw') plot(k,pacf_yw[1:naut],type='h',xaxt='n',yaxt='n',cex=0.0,cex.lab=.75,cex.axis=.75,lwd=1,xlab='',ylab='',ylim=c(-1,1))
if(method=='burg') plot(k,pacf_burg[1:naut],type='h',xaxt='n',yaxt='n',cex=0.0,cex.lab=.75,cex.axis=.75,lwd=1,xlab='',ylab='',ylim=c(-1,1))
if(method=='mle') plot(k,pacf_mle[1:naut],type='h',xaxt='n',yaxt='n',cex=0.0,cex.lab=.75,cex.axis=.75,lwd=1,xlab='',ylab='',ylim=c(-1,1))
abline(h=0)
if(limits==TRUE) {abline(h=ul,lty=2)
                    abline(h=-ul,lty=2)}
axis(side=1,cex.axis=1.4,mgp=c(3,0.6,0),tcl=-.3);
axis(side=2,las=1,cex.axis=1.4,mgp=c(3,.6,0),tcl=-.3)
if(method=='yw') mtext(side=c(1,2,1),cex=cex.labs,text=c('Lag','PACF','YW Partial Autocorrelations'),line=c(1.7,2.4,3))
if(method=='burg') mtext(side=c(1,2,1),cex=cex.labs,text=c('Lag','PACF','Burg Partial Autocorrelations'),line=c(1.7,2.4,3))
if(method=='mle') mtext(side=c(1,2,1),cex=cex.labs,text=c('Lag','PACF','MLE Partial Autocorrelations'),line=c(1.7,2.4,3))
#                }
#
if(method=='yw')out1=list(pacf_yw=pacf_yw[1:naut])
if(method=='burg')out1=list(pacf_burg=pacf_burg[1:naut])
if(method=='mle')out1=list(pacf_mle=pacf_mle[1:naut])
return(out1)                       
}

#
#

