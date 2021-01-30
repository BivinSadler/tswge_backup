plotts.dwt.wge=function (x,n.levels=4,type='S8')
{
#
#requireNamespace("waveslim", quietly = TRUE)
#
# Ben fix
#install.packages("waveslim",repos='http://cran.us.r-project.org')
#library("waveslim")
#end Ben fix

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
dj=c('d1','d2','d3','d4','d5','d6','d7','d8','d9','d10','d11','d12','d13','d14')
sj=c('s1','s2','s3','s4','s5','s6','s7','s8','s9','s10','s11','s12','s13','s14')
#
tt=1:n
if(n < 200) {plot(tt,x,type='o',xaxt='n',yaxt='n',cex=1,pch=16,cex.lab=.75,cex.axis=.75,lwd=.75,xlab='',ylab='',axes=F,xlim=c(0,n))}
if(n >= 200) {plot(tt,x,type='l',xaxt='n',yaxt='n',cex=1,pch=16,cex.lab=.75,cex.axis=.75,lwd=.75,xlab='',ylab='',axes=F,xlim=c(0,n))}
mtext(side=c(2),text='Data',cex=c(1))
#
if(type=='haar'){typewf='haar'}
if(type=='S8'){typewf='la8'}   
if(type=='D4'){typewf='d4'} 
if(type=='D6'){typewf='d6'} 
if(type=='D8'){typewf='d8'} 
x.dwt <- waveslim::dwt(x, n.levels, wf = typewf, boundary = "periodic")
x<-c(0,n)
y<-c(0,0)
ymin=0
ymax=0
nlp1=n.levels+1
for (j in 1:nlp1){
ymin=min(ymin,x.dwt[[j]])
ymax=max(ymax,x.dwt[[j]])
}
#
#
for (j in 1:n.levels)
{nj<-tt/2^j
end=n/2^j
k<-1:end
td<-k*2^j-.5*2^j+.5
plot(td,x.dwt[[j]],type='h',yaxt='n',cex=0.6,pch=15,cex.lab=2,cex.axis=2,lwd=1,xlab='',ylab='',axes=F,ylim=c(ymin,ymax),xlim=c(0,n))
lines(x,y)
mtext(side=c(2),text=dj[j],cex=c(1.))
}
plot(td,x.dwt[[n.levels+1]],type='h',yaxt='n',cex=0.6,pch=15,cex.lab=2,cex.axis=2,lwd=1,xlab='',ylab='',axes=F,ylim=c(ymin,ymax),xlim=c(0,n))
lines(x,y)
axis(side=1,cex.axis=1.5,mgp=c(3,1,0),tcl=-.3);mtext(side=c(2),text=sj[n.levels],cex=c(1.))
mtext(side=1,text='Time', cex=1.1,line=2)
if(attr(x.dwt,'wavelet') =='la8'){attr(x.dwt,'wavelet')='S8'}
if(attr(x.dwt,'wavelet') =='d4'){attr(x.dwt,'wavelet')='D4'}
if(attr(x.dwt,'wavelet') =='d6'){attr(x.dwt,'wavelet')='D6'}
if(attr(x.dwt,'wavelet') =='d8'){attr(x.dwt,'wavelet')='D8'}
return(x.dwt)
}


