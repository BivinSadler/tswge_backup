ma.pred.wge=function(x,order=3,k.ahead=1,plot=TRUE)
{


k=order
n=length(x)
nka=n+k.ahead
cat(n,k.ahead,nka,'\n')
x.sm=rep(NA,nka)
k2=k/2
k2=trunc(k2)
k2p1=k2+1
kp1=k+1
t.start=kp1
x.sm[t.start]=0
for(i in 1:k) {x.sm[t.start]=x.sm[t.start]+x[i]}
x.sm[t.start]=x.sm[t.start]/k
#cat(t.start,x.sm[t.start],'t.start,x.sm','\n')
#cat(nka,'\n')
ts1=t.start+1
ts1p1=ts1+1
np1=n+1
for(i in ts1:np1) x.sm[i]=x.sm[i-1]-x[i-1-k]/k+x[i-1]/k
#
if(k.ahead>1) {
np2=n+2
x[np1]=x.sm[np1]
for(i in np2:nka) {
x.sm[i]=0
for(j in 1:k){ 
x.sm[i]=x.sm[i]+x[i-j]/k
cat(x.sm[i],x[i-j],'\n')
}
x[i]=x.sm[i]}
                      }
if(plot==TRUE){
plot(x[1:n],type='o',xlim=c(0,nka))
points(x.sm)
points(x.sm,type='o',lwd=2)
                        }
out1 = list(x=x[1:n],pred=x.sm,order=order)
        return(out1)
}
