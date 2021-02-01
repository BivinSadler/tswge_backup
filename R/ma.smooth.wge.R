ma.smooth.wge=function(x,order=3,plot=TRUE)
{
k=order
n=length(x)
x.sm=rep(NA,n)
k2=k/2
k2=trunc(k2)
k2.chk=k2*2
if(k2.chk <k) {
k2p1=k2+1
t.start=k-k2
x.sm[t.start]=0
for(i in 1:k) {x.sm[t.start]=x.sm[t.start]+x[i]}
x.sm[t.start]=x.sm[t.start]/k
cat(t.start,x.sm[t.start],'t.start,x.sm','\n')
ts1=t.start+1
tsn=n-k2
for(i in ts1:tsn) x.sm[i]=x.sm[i-1]-x[i-k2p1]/k+x[i+k2]/k
#
}
#
#
if(k2.chk ==k) {
k2p1=k2+1
k2m1=k2-1
kp1=k+1
t.start=k-k2+1
x.sm[t.start]=x[1]/(k*2)+x[kp1]/(k*2)
for(i in 2:k) {x.sm[t.start]=x.sm[t.start]+x[i]/k}
#cat(t.start,x.sm[t.start],'t.start,x.sm','\n')
ts1=t.start+1
tsn=n-k2
for(i in ts1:tsn) x.sm[i]=x.sm[i-1]-x[i-k2p1]/(2*k)-x[i-k2]/(2*k)+x[i+k2m1]/(2*k)+x[i+k2]/(2*k)
}
if(plot==TRUE) {plot(1:n,x,type='l')
points(t.start:tsn,x.sm[t.start:tsn],type='l',lwd=2)
}
out1 = list(x=x,smooth=x.sm,order=order)
        return(out1)
}
