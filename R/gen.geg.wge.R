gen.geg.wge=function(n,u,lambda,trun=300000,vara=1,sn=0)
{
#------------------------------------------------------------------------
# Generate realization of the general linear process form
# of a Gegenbauer process based on formula (7), page 5
#                        Ferrara and Guegan 
# FORECASTING WITH K-FACTOR GEGENBAUER PROCESS: THEORY AND APPLICATIONS
# d is the lambda vector in Woodward, Gray, Zhang convention
# u is the u vector in Woodward, Gray, Zhang convention
# trun is the truncation point of the infinite GLP form to approximate
# (default for trun is 300,000)
# number is the desired length of realization
# (default for number is 200)
# mu is the mean of noise
# (default for mu is 0)
# sigma2 is the variance of noise
# (default for sigma2 is 1)
# seed is the random generation seed
# if sn=0 then produces random results, otherwise, uses seed >0
  if (sn > 0) {set.seed(sn)}
#------------------------------------------------------------------------
d=lambda 
number=n
sigma2=vara
psi=macoef.geg.wge(u,lambda,trun)
mu=0
#
# if (seed=="NULL"){ 
#                   white=rnorm(trun+number,mu,sigma2) 
#                                                       }else {set.seed(seed)
white=rnorm(trun+number,mu,sigma2)
#                       }
  realization=rep(0,number)
  for (j in 1:number)
  realization[j]=sum(psi*white[(trun+j):j])
 return(realization)
}