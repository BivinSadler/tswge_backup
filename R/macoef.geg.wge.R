macoef.geg.wge=function(u,lambda,trun=300000)
{
#------------------------------------------------------------------------
# Calculate coefficients of the general linear process form
# of a Gegenbauer process based on formula (8), page 6
#                        Ferrara and Guegan 
# FORECASTING WITH K-FACTOR GEGENBAUER PROCESS: THEORY AND APPLICATIONS
# u is the u vector in Woodward, Gray, Zhang convention
# lambda is the lambda vector in Woodward, Gray, Zhang convention
# trun is the truncation point of the infinite GLP form
# (default for trun is 300,000)
# Note: this function gives trun+1 MA coefficients for up to two GGB factors. 
#------------------------------------------------------------------------
d=lambda
 k=length(d)
 
 if (k==1)
 {
  C_0=1
  C=rep(0,trun)
  
  C[1]=2*d*u
  C[2]=2*u*((d-1)/2+1)*C[1]-(2*(d-1)/2+1)*C_0

  for (j in 3:trun)
  C[j]=2*u*((d-1)/j+1)*C[j-1]-(2*(d-1)/j+1)*C[j-2]

  psi=c(C_0,C)
 }else if (k==2)
 {
  C1_0=1
  C2_0=1
  C1=rep(0,trun)
  C2=rep(0,trun)

  C1[1]=2*d[1]*u[1]
  C1[2]=2*u[1]*((d[1]-1)/2+1)*C1[1]-(2*(d[1]-1)/2+1)*C1_0   
  for (j in 3:trun)
  C1[j]=2*u[1]*((d[1]-1)/j+1)*C1[j-1]-(2*(d[1]-1)/j+1)*C1[j-2]
 
  C2[1]=2*d[2]*u[2]
  C2[2]=2*u[2]*((d[2]-1)/2+1)*C2[2]-(2*(d[2]-1)/2+1)*C2_0   
  for (j in 3:trun)
  C2[j]=2*u[2]*((d[2]-1)/j+1)*C2[j-1]-(2*(d[2]-1)/j+1)*C2[j-2]

  psi=rep(C1_0*C2_0,trun+1)
  for (j in 2:trun+1)
  psi[j]=sum( c(C1_0,C1[1:(j-1)])*c(C2[(j-1):1],C2_0) )
 }else stop(paste("\nOnly up to two Gegenbauer factors are allowed.\n"))

 return(psi)  
}
