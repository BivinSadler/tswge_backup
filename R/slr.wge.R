slr.wge=function(x)
{
n=length(x)
time=1:n
d=lm(x~time)
#print(summary(d)$coefficients,digits=3)
z.x=x-d$coefficients[1]-d$coefficients[2]*time
#
out1=list(res=z.x,b0hat=d$coefficients[1],b1hat=d$coefficients[2],pvalue=summary(d)$coefficients[2,4],tstatistic=summary(d)$coefficients[2,3])
return(out1)
}

 
