gen.arma.wge <-
function(n,phi=0,theta=0,mu=0,vara=1,plot=TRUE,sn=0)
{
#
#  n is the realization length  (x(t), t=1, ..., n  NOTE: The generated model is based on zero mean. To generate a realization
#  with mean mu=10, for example, simply compute y(t)=x(t)+10
#  phi is a vector of AR parameters (using signs as in ATSA text)
#  theta is a vector of MA parameters (using signs as in ATSA text)  
#  vara is the white noise variance (the realization is generated with zero mean, normal white noise with variance vara)
#  plot=TRUE (default) creates a plot of the generated realization
#  NOTES:
#    (1) This function finds p and q as the length of phi and theta respectively. If either phi=0 or theta=0 (default) 
#        the values for p and q are set to zero, respectively.
#    (2) By default the white noise is zero mean, normal noise with variance vara
#    (3) This function uses a call to the base R function arima.sim which uses the same sign as ATSA for the AR parameters 
#        but opposite signs for MA parameters.  The appropriate adjustments are made here so that phi and theta should contain parameters
#        using the signs as in the ATSA text.
#        However: if you use arima.sim directly (which hs options not employed in this implementation) then you must remember that the signs 
#                 needed for the MA parameters have opposites signs as those in ATSA
# if sn=0 then produces random results, otherwise, uses seed >0
  if (sn > 0) {set.seed(sn)}
#     
#
#  This is a test
#
#  Output
#
#    Plots of the data generated
#
#    Example: For the statement test=plot3.true.atsa(n,p,phi,q,theta,vara) the following output is provided:
#             (1) a plot of the realization generated
#             (2) test$data contains the n values of the realization generated
#
#        
#
#
#Simulate Realization
#
#
ar=phi
ma=-theta
p=length(ar)
q=length(ma)
sd=sqrt(vara)
if(all(ar==0)) {ar=NA
           p=0}
if(all(ma==0)) {ma=NA
            q=0}
d=0
d1=1+d
nd=n+d
spin=2000
ngen=n+spin
#cat('n,p,d,q,ar,ma',n,p,d,q,ar,ma,'\n')
#data=arima.sim(n,model=list(order=c(p,d,q),ar=ar,ma=ma),sd=sd)
#data=arima.sim(n,model=list(ar=ar,ma=ma),sd=sd)
if((p>0) & (q>0)) {tsdata=arima.sim(ngen,model=list(order=c(p,d,q),ar=ar,ma=ma),sd=sd)}
if((p==0) & (q>0)) {tsdata=arima.sim(ngen,model=list(order=c(p,d,q),ma=ma),sd=sd)}
if((p>0) & (q==0)) {tsdata=arima.sim(ngen,model=list(order=c(p,d,q),ar=ar),sd=sd)}
if((p==0) & (q==0)) {tsdata=rnorm(ngen,mean=0,sd=sd)}
#list(order = c(1,1,0), ar = 0.7), n = 200)
y=as.numeric(tsdata)
d1=d+1
nd=n+d+1
ndspin=n+spin
xfull=rep(0,ndspin)
x=rep(0,ndspin)
for(i in d1:ndspin) {xfull[i]=y[i]}
for(ii in 1:n) {x[ii]=xfull[ii+spin]}
#
#   Plot Realization
for(ii in 1:n) {x[ii]=x[ii]+mu}
#
if(plot=='TRUE') {
cex.labs <- c(.9,.8,.9)
#
numrows <- 1
numcols <- 1
par(mfrow=c(numrows,numcols),mar=c(3.8,2.5,1,1))
t=1:n
plot(t,x[1:n],type='l',xaxt='n',yaxt='n',cex=0.5,pch=16,cex.lab=.75,cex.axis=.75,lwd=.75,xlab='',ylab='')
axis(side=1,cex.axis=.9,mgp=c(3,0.15,0),tcl=-.3);
axis(side=2,las=1,cex.axis=.9,mgp=c(3,.4,0),tcl=-.3)
mtext(side=c(1,2,1),cex=cex.labs,text=c('Time','','Realization'),line=c(1,1.1,2.1))
              }
xbar=mean(x[1:n])
varx=var(x[1:n])
return(x[1:n])     
}
