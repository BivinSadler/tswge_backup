is.sample.wge=function(data,lambda,offset)
{
#----------------------------------------------------------------------------------
# Function h.value.wge
#-----------------------------------------------------------------------------------
#-----------------------------------
# Instantaneous frequency
#------------------------------------ 
inst.freq=function(lambda,time,Gfreq)
 {
  if (lambda==0) infreq= 1/(time*(Gfreq-1)) else infreq= 1/( (time^lambda+lambda/Gfreq)^(1/lambda)-time )
  return(infreq)
 }
#----------------------------------------------
# End instantaneous frequency
#-----------------------------------------------
h.value.wge<-function(n,shift,lambda,m=n)
{
  if(lambda == 0)
    temp <- ((shift + n)/(shift + 1))^(1/(m-1))
  else temp <- (((shift + n)^lambda - (shift + 1)^lambda)/((m-1) * lambda))
  return(temp)
}
#-----------------------------------
# End h.value.wge
#------------------------------------
xlab='Time'
ylab='Frequency'
 shift=offset
dual.data=trans.to.dual.wge(data,lambda=lambda,offset=offset,h=0,plot=FALSE)
pdgm=spec.pgram(dual.data$intY,taper=0, plot=FALSE)
 per=10*log10(pdgm$spec)
 per=per-min(per)
h=h.value.wge(n=length(data),shift=shift,lambda=lambda) #n is the length of both the orig and the dual data, since they are the same length#
 t.max=length(data)
 if (lambda==0) Gfreq=h^(1/pdgm$freq) else Gfreq=pdgm$freq/h
 Gspec=per
 row.number=t.max*length(Gfreq)

 plot.matrix=matrix(0,row.number,4)  
 #[,1] is time, [,2] is Gfrequency, [,3] is calculated inst. freq, [,4] are periodogram given by pdgm, i.e. inst spec#
 plot.matrix[,1]=rep( seq(shift+1,shift+t.max),rep(length(Gfreq),t.max) )  #length(Gfreq) t=1's, then length(Gfreq.res) t=2's, until length(Gfreq) t=t.max's#
 plot.matrix[,2]=rep( Gfreq,t.max )           
 plot.matrix[,4]=rep( Gspec,t.max )

 plot.matrix[,3]=mapply(inst.freq,lambda,plot.matrix[,1],plot.matrix[,2])  
plot.matrix[,1]=plot.matrix[,1]-shift

par(mfrow=c(2,2))
plot(plot.matrix[,1],plot.matrix[,3],ylim=c(0,0.5),pch=".",cex=2,xlab=xlab,ylab=ylab,cex.lab=1.5,
      col=hcl(h=260,c=1,l=100*( 1-plot.matrix[,4]/max(plot.matrix[,4]) )) )
}



