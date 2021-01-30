plotts.mra.wge=function (x,n.levels=4,type='S8')
{
n=length(x)
dwtd=rep(0,n.levels)
   if (n/2^n.levels != trunc(n/2^n.levels)) 
        stop("Sample size is not divisible by 2^n.levels")
    if (2^n.levels > n) 
        stop("wavelet transform exceeds sample size in dwt")
numrows <- n.levels+2
numcols <- 1
timelab <- 'Lag'
valuelab <- ''
fig.width <- 5
fig.height <- 4
cex.labs <- c(1,1,1)
par(mfrow=c(numrows,numcols),mar=c(3.5,2.5,.5,1))
Sj=c('S0','S1','S2','S3','S4','S5','S6','S7','S8','S9','S10','S11','S12','S13','S14')
#
ymin=min(x)
ymax=max(x)
#
if(type=='haar'){typewf='haar'}
if(type=='S8'){typewf='la8'}   
if(type=='D4'){typewf='d4'} 
if(type=='D6'){typewf='d6'} 
if(type=='D8'){typewf='d8'} 
x.mra <- waveslim::mra(x, n.levels, wf = typewf, boundary = "periodic",method='dwt')
m=(n.levels+1)*n
sc=rep(0,m)
sm=matrix(sc,nrow=n.levels+1,ncol=n)
sm[n.levels+1,]<-x.mra[[n.levels+1]]
ymin=min(ymin,sm[n.levels+1,])
ymax=max(ymax,sm[n.levels+1,])
for(i in 1:n.levels) {
j=n.levels-i+1
sm[j,]=sm[j+1,]+x.mra[[j]]
ymin=min(ymin,sm[j,])
ymax=max(ymax,sm[j,])
}
#
# Plot data
tt=1:n
if(n < 200) {plot(tt,x,type='o',xaxt='n',yaxt='n',cex=1,pch=16,cex.lab=.75,cex.axis=.75,lwd=.75,xlab='',ylab='',axes=F,xlim=c(0,n),ylim=c(ymin,ymax))}
if(n >= 200) {plot(tt,x,type='l',xaxt='n',yaxt='n',cex=1,pch=16,cex.lab=.75,cex.axis=.75,lwd=.75,xlab='',ylab='',axes=F,xlim=c(0,n),ylim=c(ymin,ymax))}
mtext(side=c(2),text='Data',cex=c(1))
#
for(i in 1:(n.levels+1))
{plot(tt,sm[i,],type='l',yaxt='n',cex=0.6,pch=15,cex.lab=2,cex.axis=2,lwd=1,xlab='',ylab='',axes=F,ylim=c(ymin,ymax),xlim=c(0,n))
mtext(side=c(2),text=Sj[i],cex=c(1.))
}
return(x.mra)
}
