pi.weights.wge=function(phi=0,theta=0,lag.max=5) 
{
ar=phi
ma=-theta
piwt=astsa::ARMAtoAR(ar,ma,lag.max)
return(piwt)}
