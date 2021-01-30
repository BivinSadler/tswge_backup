gen.garma.wge=function(n,u,lambda,phi=0,theta=0,trun=300000,burn_in=600,vara=1,plot=TRUE,sn=0)
{
#------------------------------------------------------------------------
# Generate realization of the GLP form of a k-factor GARMA model
# based on the realization of a k-factor Gegenbauer process
# u is the u vector in Woodward, Gray, Elliott convention
# lambda is the lambda vector in Woodward, Gray, Elliott convention
# phi is the AR coefficients, INPUT 0 IF NONE
# theta is the MA coefficients, INPUT 0 IF NONE
# (default for theta is 0, i.e. no MA part)
# trun is the truncation point of the infinite GLP form to approximate
# (default for trun is 300,000)
# number is the desired length of realization
# (default for number is 200)
# burn_in is the burning-in period
# (default for burn_in is 600)
# mu is the mean of noise
# (default for mu is 0)
# vara is the variance of noise
# (default for vara is 1)
# if sn=0 then produces random results, otherwise, uses seed >0
  if (sn > 0) {set.seed(sn)}
#------------------------------------------------------------------------
number=n
mu=0 
d=lambda
sigma2=vara
if(sum(abs(d))!=0) {ggbk=length(d)} else ggbk=0
 if(sum(abs(phi))!=0) {p=length(phi)} else p=0
 if(sum(abs(theta))!=0) {q=length(theta)} else q=0
 #
 r=max(p,q)
  x=gen.geg.wge(number+burn_in+r,u,lambda,trun,sigma2,sn)

  if (p==0 & q==0) {y=x}
  else {y=rep(0,number+burn_in+r)
        for (i in (r+1):(number+burn_in+r))
        {
         y[i]=sum(phi*y[(i-1):(i-max(p,1))])+sum(c(1,-theta)*x[i:(i-max(q,1))]) 
        }
       }
#
#   Plot Realization
#
if(plot=='TRUE') {
cex.labs <- c(.9,.8,.9)
#
ymin=min(y[(1+burn_in+r):(number+burn_in+r)])
ymax=max(y[(1+burn_in+r):(number+burn_in+r)])
cat('min,max',ymin,ymax,'\n')
numrows <- 1
numcols <- 1
par(mfrow=c(numrows,numcols),mar=c(3.8,2.5,1,1))
t=1:number
plot(t,y[(1+burn_in+r):(number+burn_in+r)],type='o',xaxt='n',yaxt='n',cex=0.5,pch=16,cex.lab=.75,cex.axis=.75,lwd=.75,xlab='',ylab='',ylim=c(ymin,ymax))
axis(side=1,cex.axis=.9,mgp=c(3,0.15,0),tcl=-.3);
axis(side=2,las=1,cex.axis=.9,mgp=c(3,.4,0),tcl=-.3)
mtext(side=c(1,2,1),cex=cex.labs,text=c('Time','','Realization'),line=c(1,1.1,2.1))
              }

 return(y[(1+burn_in+r):(number+burn_in+r)])
       }

