co.wge=function(x,maxp=5)
{
maxp=5
n=length(x)
t.co=rep(0,n)
t1=1:n
# Note that b=0 in these simulations (i.e. null is true)
d=lm(x~t1)
#print(summary(d)$coefficients,digits=3)
z.x=x-d$coefficients[1]-d$coefficients[2]*t1
#cat(d$coefficients[1],d$coefficients[2],'\n')
aic.z=aic.burg.wge(z.x,p=1:maxp)
x.trans=artrans.wge(x,phi.tr=aic.z$phi,plottr=FALSE)
#AIC picks p=aic.z$p
pp=aic.z$p
#
#
p1=pp+1
#cat(p1,'p1','\n')
#cat(pp,'pp','\n')
for(tt in p1:n){
t.co[tt]=tt
sum=0
for (ii in 1:pp) {sum=sum+aic.z$phi[ii]*(tt-ii)}
t.co[tt]=tt-sum   
#cat(i,tt,aic.z$phi,t.co[tt],'i,tt.,aic.z$phi,t.co[tt]','\n')
                               }
d.co=lm(x.trans~t.co[p1:n])
#print(summary(d.co)$coefficients,digits=3)
#
#co.res=x.trans-d.co$coefficients[1]-d.co$coefficients[2]*t.co[p1:n]
#
out1=list(x=x,z.x=z.x,b0hat=d.co$coefficients[1],b1hat=d.co$coefficients[2],z.order=aic.z$p,z.phi=aic.z$phi,pvalue=summary(d.co)$coefficients[2,4],tco=summary(d.co)$coefficients[2,3])
return(out1)
}