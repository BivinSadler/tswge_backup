psi.weights.wge=function(phi=0,theta=0,lag.max=5) 
{
ar=phi
ma=-theta
psi=stats::ARMAtoMA(ar,ma,lag.max)
return(psi)
}
