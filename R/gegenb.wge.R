
#-----------------------------------------------------------------------------
# this function gegenb.wge calculates Gegenbauer polynomials (the sequence of C) 
# based on the relationship that 
#                  inf         
#                  \--
# (1-2uB+B^2)^(-d)= >  C(n;d,u)B^(n)
#                  /--
#                  n=0 
# it uses the subroutine GEGENB in the theses of Kalkomey and Chang 
# the same recursion formula can be found in Ferrara et al.
# Note: this function is the univariate version of ggb_macoeff by Ferrara and Guegan.
#-----------------------------------------------------------------------------
gegenb.wge=function(u,d,n)
{ 
 C=c( 1,2*d*u,rep(0,n-2) )
 for (i in 3:n)
 C[i]=( 2*(i-2+d)*u*C[i-1]-(i-3+2*d)*C[i-2] )/(i-1)
 return(C)
}
