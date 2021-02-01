
#
gen.aruma.wge<-function(n,phi=0,theta=0,d=0,s=0,lambda=0,vara=1, plot=TRUE,sn=0)
{
#
#  n is the realization length  (x(t), t=1, ..., n  NOTE: The generated model is based on zero #    mean. To generate a #  realization
#    with mean mu=10, for example, simply compute y(t)=x(t)+10
#  phi is a vector of AR parameters (using signs as in ATSA text)
#  theta is a vector of MA parameters (using signs as in ATSA text)  
#  d is the order of difference operator
#  s is the order of the seasonal operator
#  lambda is a vector of parameters in the nonstationary operator lambda(B) (not including the factors (1-B)^d  or (1-b^s) #  in Equation #5.8
#  vara is the white noise variance (the realization is generated with zero mean, normal white #  noise with variance vara)
#   plot=TRUE specifies to plot the generated realization
#  NOTES:
#    (1) If phi, theta, and/or lambda are known by their facrors instead of the full operator, 
#        then mult.wge should be called BEFORE calling gen.aruma.atsa
#        to create any of the vectors ph, theta, and lambda that need to be obtained
#    (2) By default the white noise is zero mean, normal noise with variance vara
#    (3) This function uses a call to the base R function arima.sim which uses the same sign 
#        as ATSA for the AR parameters 
#        but opposite signs for MA parameters.  The appropriate adjustments are made here so 
#        that phi and theta should contain parameters
#        using the signs as in the ATSA text.
#        However: if you use arima.sim directly then you must remember that the signs needed 
#        for the MA parameters have opposites signs as those in ATSA
#    (4) gen.aruma.atsa is a generalization of the functions gen.arma.atsa and 
#        gen.arima.atsa, and it can be the only generating function that you use if you so 
#        choose by correctly specifying parameters.
# if sn=0 then produces random results, otherwise, uses seed >0
  if (sn > 0) {set.seed(sn)}
#
#
#
#
#  Output
#
#    Plots of the data generated
#        
cex.labs <- c(.9,.8,.9)
#
numrows <- 1
numcols <- 1
par(mfrow=c(numrows,numcols),mar=c(3.8,2.5,1,1))
#
#
#
#Simulate Realization
#
#
#cat('ar, ma',ar,ma,'\n')
ar=phi
ma=-theta
p=length(ar)
q=length(ma)
dlam=length(lambda)
sd=sqrt(vara)
if(all(ar==0)) {ar=NA
           p=0}
if(all(ma==0)) {ma=NA
            q=0}
if(all(lambda==0)) {lambda=NA
            dlam=0}
sd=sqrt(vara)
dlams=dlam+s
lambdas=rep(0,100)
lambdas
seas=rep(0,100)
if(s>0) seas[s]=1
spin=100
ngen=n+dlams+spin
if((p>0) & (q>0)) {tsdata=arima.sim(ngen,model=list(order=c(p,d,q),ar=ar,ma=ma),sd=sd)}
if((p==0) & (q>0)) {tsdata=arima.sim(ngen,model=list(order=c(p,d,q),ma=ma),sd=sd)}
if((p>0) & (q==0)) {tsdata=arima.sim(ngen,model=list(order=c(p,d,q),ar=ar),sd=sd)}
if((p==0) & (q==0)) {tsdata=arima.sim(ngen,model=list(order=c(0,d,0)),sd=sd)}
#
tsdata
#
#  Compute the inverse of the nonstationary operator  (similar to diffinv in R)
#
#  
y=as.numeric(tsdata)
if ((dlam>0) & (s>0)) {temp=mult.wge(fac1=lambda,fac2=seas)
                       lambdas=temp$model.coef}
if ((dlam>0) & (s==0)) {lambdas=lambda}
if ((dlam==0) & (s>0)) {lambdas=seas}
if ((dlam==0) & (s==0)) {lambdas=0}

lambdas
d1=d+dlams+1
nd=n+d+dlams+1
ndspin=nd+spin-1
xfull=rep(0,ndspin)
x=rep(0,ndspin)
if(dlams == 0) {
for (i in d1:ndspin) {xfull[i]=y[i]}
                }
if(dlams > 0) {
for (i in d1:ndspin) {
xfull[i]=y[i]
for (j in 1:dlams) {
xfull[i]=xfull[i]+lambdas[j]*xfull[i-j]
                    }
                       }
                 }
for (ii in 1:n) {
x[ii]=xfull[ii+spin+d1-1]
}


#
#   Plot Realization
#

if(plot=='TRUE') {t=1:n
x=as.numeric(x)
plot(t,x[1:n],type='o',xaxt='n',yaxt='n',cex=0.5,pch=16,cex.lab=.75,cex.axis=.75,lwd=.75,xlab='',ylab='')
axis(side=1,cex.axis=.9,mgp=c(3,0.15,0),tcl=-.3);
axis(side=2,las=1,cex.axis=.9,mgp=c(3,.4,
0),tcl=-.3)
mtext(side=c(1,2,1),cex=cex.labs,text=c('Time','','Realization'),line=c(1,1.1,2.1))
}

return(x[1:n])         
}

