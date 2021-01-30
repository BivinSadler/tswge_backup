
#
#
#
#   For a given ARMA(p,q) model,  this function generates a realization, calculates the true autocorrelations 
#   and spectral density (dB) and plots the three graphs (similar to Figure 3.10 in the ATSA text
#
plotts.true.wge<-function(n=100,phi=0,theta=0,lag.max=25,vara=1)
{
#
#  n is the realization length  (x(t), t=1, ..., n. (Default n=100)  NOTE: The generated model is based on zero mean. To generate a realization
#  with mean mu=10, for example, simply compute y(t)=x(t)+10
#  phi is a vector of AR parameters (using signs as in ATSA text)  NOTE:  If p=0 then phi must be 0 also
#  theta is a vector of MA parameters (using signs as in ATSA text)  
#  lag.max is the maximum lag at which the autocorrelations and autocovariances will be calculated
#  vara is the white noise variance
# 
#  NOTEs:
#    (1) This function finds p and q as the length of phi and theta respectively. If either phi=0 or theta=0 (default) 
#        the values for p and q are set to zero, respectively.
#    (2) By default the white noise is zero mean, normal noise with variance vara
#    (3) max(p,q+1)<=25
#    (4) This function uses a call to the base R function arima.sim which uses the same sign as ATSA for the AR parameters 
#        but opposite signs for MA parameters.  The appropriate adjustments are made here so that phi and theta should contain parameters
#        using the signs as in the ATSA text.
#        However: if you use arima.sim directly (which hs options not employed in this implementation) then you must remember that the signs 
#                 needed for the MA parameters have opposites signs as those in ATSA
#
#
#  Output
#
#    Plots of the data generated
#
#    Example: For the statement test=plot3tr.wge(n,p,phi,q,theta,lag.max,vara) the following output is provided:
#             test$data contains the n values of the realization generated
#             test$aut1 contains the true autocorrelations (lags 0, 1, ..., lag.max)
#             test$g contains the true autocovariances (lags 0, 1, ..., lag.max)  Note: test$g[1] contains the process variance
#             test$gvar also provides the process variance, i.e. gvar=g[1]
#             test#spec contain the 251 spectral density values associated with f=0, 1/500, 2/500, ..., 250/500              
p=length(phi)
q=length(theta)
sd=sqrt(vara)
if(all(phi==0)) {ar=NA
           p=0}
if(all(theta==0)) {ma=NA
            q=0}
nm1=n-1
if(lag.max >nm1) {lag.max=nm1}
ip=max(p,q+1)
ipm1=ip-1
ipp1=ip+1
lag.maxp1=lag.max+1
g=rep(0,lag.maxp1)
d=rep(0,lag.maxp1)
aut1=rep(0,lag.maxp1)
aut1
spec=rep(0,251)
a=matrix(rep(0,676), ncol = 26)
papp=ip-p
#
#
#   Set up for plots
#
#
layout(mat=matrix(  c(0,1,1,0,
                      2,2,3,3)   ,nrow=2,byrow=TRUE))
par(mar=c(5,2.8,1,1.5))
cex.labs <- c(.9,.8,.9)
#
#
#Simulate Realization
#
#
#cat('ar, ma',ar,ma,'\n')
ar=phi
ma=-theta
if((p>0) & (q>0)) {data=arima.sim(n,model=list(ar=ar,ma=ma))}
if((p==0) & (q>0)) {data=arima.sim(n,model=list(ma=ma))}
if((p>0) & (q==0)) {data=arima.sim(n,model=list(ar=ar))}
if((p==0) & (q==0)) {data=arima.sim(n)}
t1=1:length(data)
plot(t1,data,type='l',xaxt='n',yaxt='n',cex=0.75,pch=15,cex.lab=.75,cex.axis=.75,lwd=.75,xlab='',ylab='')
axis(side=1,cex.axis=.85,mgp=c(3,0.15,0),tcl=-.3);
axis(side=2,las=1,cex.axis=.85,mgp=c(3,.4,0),tcl=-.3)
mtext(side=c(1,2,1),cex=cex.labs,text=c('Time','','(a) Realization'),line=c(1.25,1.4,2.5),font=c(1,1,1))
#
#
#  Calculate True Autocorrelations
#
#
PHI<-phi
  while(papp>0) {zp<-rep(0,papp)
                 PHI<-append(phi,zp)
                 papp=-1
                 }
PHI
qapp=ip-q
qapp
THETA<-theta
  while(qapp>0) {zq=rep(0,qapp)
                THETA=append(theta,zq)
                qapp=-1
                 }
THETA
#
if(ip<=1) {d1=(PHI[1]-THETA[1])*vara
           dn=1-PHI[1]^2           
           g[1]=(vara-THETA[1]*d1-THETA[1]*vara*PHI[1])/dn
           }
#
#
#
if(ip >1) {ipm1=ip-1
           ipp1=ip+1
           d[1]=vara
# do 10
           for (i in 2:ip) {
                            im1=i-1
# do 20
                            for (j in 1:im1) {
                                              d[i]=d[i]+PHI[j]*d[i-j]
                                              }
                             d[i]=d[i]-THETA[i-1]*vara
                            }
# do 40
            for (i in 1:ipm1) {
                               g[i]=0
                               k=ip-i
# do 40
                               for (j in 1:k) {
                                               g[i]=g[i]-THETA[j+i-1]*d[j+1]
                                               }                 
                                }
              g[ip]=0
              g[1]=g[1]+d[1]
# do 50
              for (i in 2:ip) {
                               g[i]=g[i]-THETA[i-1]*d[1]
                               }
# do 60
              for (i in 1:ip) {
                               a[i,i]=1
                               }
# do 80
              for (i in 1:ipm1) {
                                 for (j in 1:ip) {
                                                  ii=abs(i-j+1)+1
                                                  a[j,ii]=a[j,ii]-PHI[i]
                                                  }
                                 }
# do 90
              for (i in 2:ip) {
                               ii=abs(ip-i+1)+1
                                a[i,ii]=a[i,ii]-PHI[ip]
                               }
# do 100
               for (i in 1:ip) {
                               ii=ip-i+1
                               a[1,i]=a[1,i]-PHI[ip]*PHI[ii]
                               }
# do 120
                for(k in 2:ip) {
                                m=k-1
                                
# do 120
                                 for (i in k:ip) {
# do 130                                        
                                                 for (j in k:ip) {
                                                                 a[i,j]=a[i,j]-a[i,m]*a[m,j]/a[m,m]
                                                                 }
                                                  g[i]=g[i]-a[i,m]*g[m]/a[m,m]
                                                  }
                                 }
                  g[ip]=g[ip]/a[ip,ip]
# do 140
                  for (i in 1:ipm1) {
                                     m=ip-i
                                     ctc=0
# do 150
                                     for (j in 1:i) {
                                                     k=ip-j+1
                                                     ctc=ctc+a[m,k]*g[k]
                                                     }                                      
                                      g[m]=(g[m]-ctc)/a[m,m]
                                      }
           }


#
#
#use difference equation to complete the list of autocorrelations
#
#
gvar=g[1]
lag.maxp1=lag.max+1
for (i in 1:lag.maxp1) {aut1[i]=g[i]/gvar}



if(p > 0) {

           for (k in ipp1:lag.maxp1) {
                                 g[k]=0
                                 for (j in 1:p) {
                                                 g[k]=g[k]+PHI[j]*g[k-j]                                                }
                                                 }
           }
gvar=g[1]
for (i in 1:lag.maxp1) {aut1[i]=g[i]/gvar
                                      }
#
#
# Calculate Spectral Density
#
#
#
I=sqrt(as.complex(-1))
Pi=3.14159265359
for (fi in 1:251) {
                  num=1
                  den=1
                  f=(fi-1)/500
                  if(q>0) {
                           for (k in 1:q) {
                            num=num-theta[k]*exp(-2*Pi*I*k*f)
                                           }
                  if(q==0) {num=1}
                           }
                  if (p>0) {                  
                            for (k in 1:p) {
                            den=den-phi[k]*exp(-2*Pi*I*k*f)
                                            }
                  if(p==0) {den=1}
                            }
                   spec[fi]=10*log10((vara*abs(num)^2)/(gvar*abs(den^2)))
#cat('f,num,den,spec(f)',f,num, den, spec[fi],'\n')
                 }
#
# plot true autocorrelations
#
k=0:lag.max
plot(k,aut1,type='h',xaxt='n',yaxt='n',cex=0.4,cex.lab=.75,cex.axis=.75,lwd=.75,xlab='',ylab='',ylim=c(-1,1))
abline(h=0)
axis(side=1,cex.axis=.9,mgp=c(3,0.15,0),tcl=-.3);
axis(side=2,las=1,cex.axis=.9,mgp=c(3,.4,0),tcl=-.3)
mtext(side=c(1,2,1),cex=cex.labs,text=c('Lag','','(b) True Autocorrelations'),line=c(1,1.1,2.1))
#
# plot true spectral density
#
fii=1:251
f=(fii-1)/500
plot(f,spec,type='l',xaxt='n',yaxt='n',cex=0.4,cex.lab=.75,cex.axis=.75,lwd=.75,xlab='',ylab='')
axis(side=1,cex.axis=.9,mgp=c(3,0.15,0),tcl=-.3);
axis(side=2,las=1,cex.axis=.9,mgp=c(3,.4,0),tcl=-.3)
mtext(side=c(1,2,1),cex=cex.labs,text=c('Frequency','dB','(c) True Spectral Density'),line=c(1,1.1,2.1))
#
out1=list(data=data,aut1=aut1,acv=g,spec=spec)
return(out1)
}

