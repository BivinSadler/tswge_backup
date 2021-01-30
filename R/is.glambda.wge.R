
#---------------------------------------------------------------------------------------------------------
# glambda IS plot
# based on AR model fit to dual data
#---------------------------------------------------------------------------------------------------------
is.glambda.wge=function(n,phi=0,sigma2=1,lambda,offset)
{
shift=offset
xlab='Time'
ylab='Frequency'
theta=0
freq=seq(0.001,.500,.001)
t.max=n
spec=c() 
u=1
d=0
#--------------------------------------------------------------------------
# Calculate the theoretical spectrum of AR model at frequency f
#--------------------------------------------------------------------------
spec.ar.f.wge=function(f,phi,sigma2=1){
if(sum(abs(phi))!=0) {p=length(phi)} else p=0
# phi{exp(i*2*pi*f)}
phi_func=sum(c(1,-phi)*exp(2*pi*1i*f)^c(0:p))
phi_norm=(Mod(phi_func))^2
spect=sigma2*(1/phi_norm)
return(spect)
}
# ------------------------------------------------------------------------------
#End spec.ar.f.wge
#-------------------------------------------------------------------------------
#
#----------------------------------------------------------------------------------
# Function h.value.wge
#-----------------------------------------------------------------------------------
h.value.wge<-function(n,shift,lambda,m=n)
{
  if(lambda == 0)
    temp <- ((shift + n)/(shift + 1))^(1/(m-1))
  else temp <- (((shift + n)^lambda - (shift + 1)^lambda)/((m-1) * lambda))
  return(temp)
}
#-------------------------------------------------------------------------------------------------
# End h.value.wge
#-------------------------------------------------------------------------------------------------
#
# Function inst.freq
#-------------------------------------------------------------------------------------------------
inst.freq=function(lambda,time,Gfreq)
 {
  if (lambda==0) infreq= 1/(time*(Gfreq-1)) else infreq= 1/( (time^lambda+lambda/Gfreq)^(1/lambda)-time )
  return(infreq)
 }
#--------------------------------------------------------------------------------------
# End inst.freq
#-------------------------------------------------------------------------------------
#
#--------------------------------------------------------------
#  reinterpol.wge
#--------------------------------------------------------------
reinterpol.wge<-function(data,start,h,lambda,int.start,length,IntMethod="Linear")
{
  if(lambda == 0) {
    m <- log(start)/log(h)
    ht <- h^seq(m, (m + length(data)-1))
  }
  else {
    m <- (start^lambda - 1)/(h * lambda)
    kk <- seq(m, (m + length(data)-1))
    ht <- (kk * h * lambda + 1)^(1/lambda)
  }
  tt <- seq(int.start, (int.start + length-1))
  if(IntMethod=="Linear") intData <- approx(ht, data, xout = tt, rule = 3)$y
   return(intData)
}
#--------------------------------------------------------------
#  end reinterpol.wge
#--------------------------------------------------------------
#
#
for (i in 1:length(freq)) {spec[i]=spec.ar.f.wge(freq[i],phi,sigma2)}
spec=10*log10(spec) 
spec=spec-min(spec)
h=h.value.wge(n=t.max,shift=offset,lambda=lambda)
 if (lambda==0) Gfreq=h^(1/freq) else Gfreq=freq/h
 Gspec=spec
 row.number=t.max*length(Gfreq)
#
 plot.matrix=matrix(0,row.number,4)  
 #[,1] is time, [,2] is Gfrequency, [,3] is calculated inst. freq, [,4] are theoretic dual model spec, i.e. inst spec#
 plot.matrix[,1]=rep( seq(shift+1,shift+t.max),rep(length(Gfreq),t.max) )  #length(Gfreq) t=1's, then length(Gfreq.res) t=2's, until length(Gfreq) t=t.max's#
 plot.matrix[,2]=rep( Gfreq,t.max )           
 plot.matrix[,4]=rep( Gspec,t.max )
#
 plot.matrix[,3]=mapply(inst.freq,lambda,plot.matrix[,1],plot.matrix[,2])  
plot.matrix[,1]=plot.matrix[,1]-shift
par(mfrow=c(2,2))
 plot(plot.matrix[,1],plot.matrix[,3],ylim=c(0,0.5),pch=".",cex=2,xlab=xlab,ylab=ylab,cex.lab=1.5,cex.main=1.5,
      col=hcl(h=260,c=1,l=100*( 1-plot.matrix[,4]/max(plot.matrix[,4]) )) )
}

