#
#
#
#
true.arma.spec.wge<-function(phi=0,theta=0,vara=1,plot=TRUE)
{
#
#  phi is a vector of AR parameters (using signs as in ATSA text)
#  theta is a vector of MA parameters (using signs as in ATSA text)  
#   plot is a logical variable: TRUE=plot spectral density, FALSE=no plot}
# 
#  NOTEs:
#    (1) By default the white noise is zero mean, normal noise with variance vara
#    (2) max(p,q+1)<=25
#    (3) This function uses a call to the base R function arima.sim which uses the same sign as ATSA for the AR parameters 
#        but opposite signs for MA parameters.  The appropriate adjustments are made here so that phi and theta should contain parameters
#        using the signs as in the ATSA text.
#        However: if you use arima.sim directly (which hs options not employed in this implementation) then you must remember that the signs 
#                 needed for the MA parameters have opposites signs as those in ATSA
#
#
#  Output
#
#    Plots of the spectral density generated
#
#    Example: For the statement test=true.arma.spec.wge(phi,theta,vara=1) the following output is provided:
#             test$spec contain the 251 spectral density values associated with f=0, 1/500, 2/500, ..., 250/500              
#
p=length(phi)
q=length(theta)
if(all(phi==0)) {ar=NA
           p=0}
if(all(theta==0)) {ma=NA
            q=0}
ip=max(p,q+1)
ipm1=ip-1
ipp1=ip+1
maxk=ip
maxkp1=maxk+1
g=rep(0,maxkp1)
d=rep(0,maxkp1)
aut1=rep(0,maxkp1)
#vara=1
spec=rep(0,251)
a=matrix(rep(0,676), ncol = 26)
papp=ip-p
#
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
maxkp1=maxk+1
for (i in 1:maxkp1) {aut1[i]=g[i]/gvar}



if(p > 0) {

           for (k in ipp1:maxkp1) {
                                 g[k]=0
                                 for (j in 1:p) {
                                                 g[k]=g[k]+PHI[j]*g[k-j]                                                }
                                                 }
           }
gvar=g[1]
for (i in 1:maxkp1) {aut1[i]=g[i]/gvar
                                      }
#
#
#
#
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
if(plot=='TRUE') {
cex.labs <- c(.9,.8,.9)
#numrows <- 1
#numcols <- 1
#par(mfrow=c(numrows,numcols),mar=c(3.8,2.5,1,1))
#
# plot true spectral density
#
fii=1:251
f=(fii-1)/500
plot(f,spec,type='l',xaxt='n',yaxt='n',cex=0.4,cex.lab=.75,cex.axis=.75,lwd=.75,xlab='',ylab='')
axis(side=1,cex.axis=.9,mgp=c(3,0.15,0),tcl=-.3);
axis(side=2,las=1,cex.axis=.9,mgp=c(3,.4,0),tcl=-.3)
mtext(side=c(1,2,1),cex=cex.labs,text=c('Frequency','dB','True Spectral Density'),line=c(1,1.1,2.1))
#
}
out1=list(f=f,spec=spec)
return(out1)
}
