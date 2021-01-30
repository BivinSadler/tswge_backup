psi.weights.wge=function(phi=0,theta=0,lag.max=0) 
{
ar=phi
ma=-theta
psi=ARMAtoMA(ar,ma,lag.max)
return(psi)
}
