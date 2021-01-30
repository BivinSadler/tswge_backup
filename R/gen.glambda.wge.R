#############################
#  generate G-lambda sample
#############################
gen.glambda.wge=function(n,lambda,phi=0,offset=20,vara=1,plot=TRUE,sn=0)
{
#------------------------------------------------
#function h.value.wge
#------------------------------------------------
h.value.wge<-function(n,shift,lambda,m=n)
{
  if(lambda == 0)
    temp <- ((shift + n)/(shift + 1))^(1/(m-1))
  else temp <- (((shift + n)^lambda - (shift + 1)^lambda)/((m-1) * lambda))
  return(temp)
}
#---------------------------------------------------
# end function h.value.wge
#---------------------------------------------------

IntMethod='Linear'
ar=phi
ma=0
shift=offset
var=vara
mean=0
h=0
#  if (sn > 0) {set.seed(sn)}
  eq.x <- seq((shift + 1), (shift + n))
  if(lambda == 0) {
    
    if(h <= 1) {
      h <- h.value.wge(n=n,shift=shift,lambda=0)
      k<-log(shift+1)/log(h)
      eu.x <- h^seq(k,(k+n-1))
      if(sn > 0) set.seed(sn)
      yt <- (arima.sim(n = n+ 600, n.start = 600, model = list(ar = ar,ma = ma))[
        -c(1:600)]) * sqrt(var)+mean
    }
    else {
      n.h <- floor(log(n + shift)/log(h))+1
      eu.x <- h^seq(0, n.h)
      if(sn > 0)
        set.seed(sn)
      yt <- (arima.sim(n = n.h + 601, n.start = 600, model = list(ar = ar,ma = ma))[ 
        -c(1:600)]) * sqrt(var)+mean
    }
  } else {
    if(h <= 0){
      h<-h.value.wge(n=n,shift=shift,lambda=lambda)
      k<-((shift+1)^lambda-1)/(h*lambda)
      eu.x <- (h*lambda*seq(k, (k+n-1)) + 1)^(1/lambda)
      if(sn > 0) set.seed(sn)
      yt <- (arima.sim(n = n+ 600, n.start = 600, model = list(ar = ar,ma = ma))[ 
        -c(1:600)]) * sqrt(var)+mean
    } else {
      n.h <- floor(((n + shift)^lambda - 1)/(h*lambda))+1
      eu.x <- (h* lambda * seq(0, n.h) + 1)^(1/lambda)
      if(sn > 0)  set.seed(sn)
      yt <- (arima.sim(n = n.h + 601, n.start = 600, model = list(ar = ar))[ 
        -c(1:600)]) * sqrt(var)+mean
    }
  }
  #equal space & interpolate
  
  if(IntMethod == "Linear")       
    into <- approx(eu.x, yt, xout = eq.x, rule = 3)$y
  else stop(paste("\nWe don't have other methods yet\nPlease selec\nt oneof the previous two\n"))
  
#return(list(TVFdata = into, oriX = eu.x[eu.x > shift ], oriY =yt[eu.x >shift],h=h))
tvfdata=into
#
#   Plot Realization
#
if(plot=='TRUE') {
cex.labs <- c(.9,.8,.9)
#
numrows <- 1
numcols <- 1
par(mfrow=c(numrows,numcols),mar=c(3.8,2.5,1,1))
t=1:n
plot(t,tvfdata,type='o',xaxt='n',yaxt='n',cex=0.5,pch=16,cex.lab=.75,cex.axis=.75,lwd=.75,xlab='',ylab='')
axis(side=1,cex.axis=.9,mgp=c(3,0.15,0),tcl=-.3);
axis(side=2,las=1,cex.axis=.9,mgp=c(3,.4,0),tcl=-.3)
mtext(side=c(1,2,1),cex=cex.labs,text=c('Time','','Realization'),line=c(1,1.1,2.1))
              }
return(tvfdata)
}
