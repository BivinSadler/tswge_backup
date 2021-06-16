mult.wge=function(fac1=0,fac2=0,fac3=0,fac4=0,fac5=0,fac6=0) 
{
fac1=-fac1
fac2=-fac2
fac3=-fac3
fac4=-fac4
fac5=-fac5
fac6=-fac6
bfac1=c(1,fac1)
bfac2=c(1,fac2)
bfac3=c(1,fac3)
bfac4=c(1,fac4)
bfac5=c(1,fac5)
bfac6=c(1,fac6)
##make into poynomials##   pfaci
pfac1=PolynomF::polynom(bfac1)
pfac2=PolynomF::polynom(bfac2)
pfac3=PolynomF::polynom(bfac3)
pfac4=PolynomF::polynom(bfac4)
pfac5=PolynomF::polynom(bfac5)
pfac6=PolynomF::polynom(bfac6)
#
#  Multiply polynomials
#
mfac=pfac1*pfac2*pfac3*pfac4*pfac5*pfac6
#
#
#extract coefficients
#
tmpcoef=coef(mfac)
#
#
#
#remove first element
#
# if all factors are 0 then return 0...this can happen if fore.arima.wge is called with d = 0
if(length(tmpcoef) > 1)
{
model.coef=tmpcoef[-1]
}
else
{
  model.coef = 0
}
#
#
#correct signs
model.coef=-model.coef
#
#
out1=list(char.poly=mfac,model.coef=model.coef)
return(out1)  
}




