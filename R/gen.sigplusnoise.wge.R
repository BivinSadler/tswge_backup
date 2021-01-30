#  Program to generate signal + noise realizations with a maximum of two cosine terms
#  
#
#
#
#
#
gen.sigplusnoise.wge<-function(n,b0=0,b1=0,coef=c(0,0),freq=c(0,0),psi=c(0,0),phi=0,vara=1,plot='TRUE',sn=0)
{
#
# x(t)=b0+b1*t+coef[1]*cos(2*pi*freq[1]*t+psi[1])+coef[2]*cos(2*pi*freq[2]*t+psi[2])+z(t)
#
#   Input:
#   n=length of realization to be generated  (t in the above formula spacify the integers from 1 to n)
#   coef is a 2-component vector specifying the coefficients (if only one cosine term is desired define coef[2]=0)
#   freq is a 2-component vector specifying the frequency components (0 to .5)
#   psi is a 2-component vector specifying the phase shift (0 to 2pi)
#   z(t) is a zero mean realization from the AR model phi(B)z(t)=a(t)   where a(t) is N(0,vara) white noise
# 
#   plot is a logical variable. 'TRUE' produces a plot of the generated signal plus noise data set.
#
# if sn=0 then produces random results, otherwise, uses seed >0
  if (sn > 0) {set.seed(sn)}#
#
#Output:
#

#
#
#
std=sqrt(vara)
t=1:n
zt=gen.arma.wge(n=n,phi=phi,theta=0,vara=vara,sn=sn)
x=b0+b1*t+coef[1]*cos(2*pi*freq[1]*t+psi[1])+coef[2]*cos(2*pi*freq[2]*t+psi[2])+zt
numrows <- 1
numcols <- 1
fig.width <- 5.5
fig.height <- 4.5
cex.labs <- c(.8,.7,.9)
par(mfrow=c(numrows,numcols),mar=c(3.8,2.5,1,1))
if(plot =='TRUE'){ plot(t,x,type='l')}
return(x)
}
