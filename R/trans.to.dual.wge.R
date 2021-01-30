trans.to.dual.wge=function(x,lambda,offset=60,h=0,plot=TRUE)
{
#------------------------------------------------------------
# function h.value.wge
#------------------------------------------------------------
h.value.wge<-function(n,shift,lambda,m=n)
{
  if(lambda == 0)
    temp <- ((shift + n)/(shift + 1))^(1/(m-1))
  else temp <- (((shift + n)^lambda - (shift + 1)^lambda)/((m-1) * lambda))
  return(temp)
}
#------------------------------------------------------------
# end function h.value.wge
#------------------------------------------------------------
#
data=x
n=length(data)
shift=offset 
IntMethod="Linear"
if((h<=1&&lambda==0)||(h<=0&&lambda!=0)){h=h.value.wge(n, shift, lambda,m=n)}
 tt=seq((shift+1),(shift+n))
 if(lambda==0){
  	         k=log(shift + 1)/log(h) 
		   m=log(shift + n)/log(h)
		   kk=seq(k, m)
		   ht2=h^kk
	        }
 else if(lambda==1){ht2 <- seq((shift + 1),(shift + n))}
 else {k=((1+shift)^lambda-1)/(h*lambda) 
	 m=((n+shift)^lambda-1)/(h*lambda)
	 kk=seq(k, m)
	 ht2=(kk*h*lambda+1)^(1/lambda)
	}
	
 if(IntMethod=="Linear"){
		             intData=approx(tt,data,xout=ht2,rule=2)
		             intY=intData$y[round(ht2,2)>=(shift+1)]
	                  }
 else stop(paste("\nWe don't have other methods yet\nPlease select one of the previous three\n"))

#   Plot Realization
#
if(plot=='TRUE') {
cex.labs <- c(.9,.8,.9)
#
numrows <- 1
numcols <- 1
par(mfrow=c(numrows,numcols),mar=c(3.8,2.5,1,1))
t=1:n
plot(t,intY,type='o',xaxt='n',yaxt='n',cex=0.5,pch=16,cex.lab=.75,cex.axis=.75,lwd=.75,xlab='',ylab='')
axis(side=1,cex.axis=.9,mgp=c(3,0.15,0),tcl=-.3);
axis(side=2,las=1,cex.axis=.9,mgp=c(3,.4,0),tcl=-.3)
mtext(side=c(1,2,1),cex=cex.labs,text=c('Time','','Realization'),line=c(1,1.1,2.1))
              }
list(intX=ht2[round(ht2,2)>=(shift+1)],intY=intY,h=h)
}
