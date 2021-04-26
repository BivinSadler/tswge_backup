wbg.boot.wge=function(x,nb=399,alpha=.05,pvalue=TRUE,sn=0){
# Finds observed significance level using co.boot
  if (sn > 0) {set.seed(sn)}
n=length(x)
xb.t.co=rep(0,399)
xb.t.co.abs=rep(0,399)
w=co.wge(x,maxp=5)
x.t.co=w$tco
x.aic=aic.burg.wge(x,p=1:5)
for(ii in 1:nb){ xb=gen.arma.wge(n,phi=x.aic$phi,plot=FALSE)
wb=co.wge(xb,maxp=5)
xb.t.co[ii]=wb$tco
#cat(ii,xb.t.co[ii],'ii,xtbco[ii]','\n')
xb.t.co.abs[ii]=abs(xb.t.co[ii])  }
y=sort(xb.t.co.abs)
#crit=(1-alpha)*(nb+1)
#k=0
#if(abs(x.t.co)>y[crit]) {k=k+1}
f=0
for (tt in 1:nb) {if(y[tt]>=abs(x.t.co)) f=f+1}
pv=f/nb
out1=list(p=x.aic$p,phi=x.aic$phi,pv=pv)
return(out1)
}

