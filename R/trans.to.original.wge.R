trans.to.original.wge=function(xd,lambda,offset,h,plot=TRUE)
{
IntMethod="Linear" 
data=xd
length=length(xd)
n=length(xd)
start=offset+1
int.start=offset
if(lambda==0){
		   m=log(start)/log(h)
		   ht=h^seq(m,(m+length(data)-1))
	        }
 else{
	m=(start^lambda-1)/(h*lambda)
	kk=seq(m,(m+length(data)-1))
	ht=(kk*h*lambda+1)^(1/lambda)
     }
 tt=seq(int.start,(int.start+length-1))
	
 if(IntMethod=="Linear"){intData=approx(ht,data,xout=tt,rule=2)$y}  #ht is the timescale of all the data points#
                                                                    #data is the full length data#
                                                                    #xout or tt can be only a few time points#
                                                                    #then approx only interpolates for those time points#
                                                                    #this is what happens in the forecast program# 
 else stop(paste("\nWe don't have other methods yet\nPlease select one of the previous two\n"))
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
plot(t,intData,type='o',xaxt='n',yaxt='n',cex=0.5,pch=16,cex.lab=.75,cex.axis=.75,lwd=.75,xlab='',ylab='')
axis(side=1,cex.axis=.9,mgp=c(3,0.15,0),tcl=-.3);
axis(side=2,las=1,cex.axis=.9,mgp=c(3,.4,0),tcl=-.3)
mtext(side=c(1,2,1),cex=cex.labs,text=c('Time','','Realization'),line=c(1,1.1,2.1))
              }
 return(intData)
}
