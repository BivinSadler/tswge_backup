#----------------------------------------------------------------------
# Forecast based on grid-search-estimated GARMA model
#----------------------------------------------------------------------
fore.garma.wge=function(x,u,lambda,phi=0,theta=0,n.ahead=10,lastn=TRUE,plot=TRUE)
{
# the number of Gegenbauer factors in the proposed model
k=length(u)
d=lambda
forecast.step=n.ahead
nback=500
if(lastn==TRUE){forecast.step=-forecast.step}
if(sum(abs(phi))!=0) {p=length(phi)} else p=0
if(sum(abs(theta))!=0) {q=length(theta)} else q=0


# the number of observations in the input realization 
narma=length(x)

# the sequence of backcasted observations (for now all are 0) 
# and input realization (centered) 
nc=narma+nback
z=rep(0,nc)
z[-(1:nback)]=x-mean(x)

# fit a high order AR model for realization to backcast
# use ar.yw since, after specifying AR order, Yule-Walker gives the optimal 
# AR coefficient estimates in AR model case.
AR.p=20
AR.phi=ar.yw(z[-(1:nback)],order.max=AR.p,aic=FALSE,demean=FALSE)$ar

# use this high order AR model to backcast X(0),...X(-nback+1) 
# that is, correspondingly, z(nback),...z(1) 
for (i in nback:1) {z[i]=sum(AR.phi*z[(i+1):(i+AR.p)])}

if (k==1)
{           # transform z sequence to remove Gegenbauer factors
            # and derive a new sequence of w, which can be fit into ARMA model
            C1=gegenb.wge(u,-d,nc)
            w=c()
            for (i in 1:narma) {w[i]=sum(C1[1:(nback+i-1)]*z[(nback+i):2])}
            w=w-mean(w)
            w.forecast=w
  # if forecast.step>0, forecast forecast.step steps in the future
  # if forecast.step<0, update the last forecast.step steps with fitted AR model
  # meanwhile fit x with an AR model and compare the forecast made by AR and GARMA models

  x.ARest=ar.yw(x,aic=TRUE,order.max=8)
if(x.ARest$order>0) {ma.constant=1
for(i in 1:x.ARest$order){ma.constant=ma.constant-x.ARest$ar[i]}
ma.constant=ma.constant*mean(x)
}
  x.ARfore=x
  if (forecast.step<0) 
{
 for (i in (narma+forecast.step+1):narma) {w.forecast[i]=sum(w.forecast[(i-1):(i-p)]*phi)
                                             x.ARfore[i]=ma.constant+sum(x.ARfore[(i-1):(i-x.ARest$order)]*x.ARest$ar)
                                            }
  }
  if (forecast.step>0) 
  {
   w.forecast=c(w.forecast,predict(arima(w.forecast,order=c(p,0,q)),n.ahead=forecast.step)$pred)
   x.ARfore=c(x.ARfore,predict(x.ARest,n.ahead=forecast.step)$pred)
#cat('w.forecast,x,ARfore',w.forecast,x.ARfore,'\n')
  }
  # so now the w.forecast either still has narma points (if last forecast.step steps updated)
  # or it has narma+forecast.step points (if forecast.step steps in future are forecasted)
  narma2=length(w.forecast)
  nc2=nback+narma2
  z2=rep(0,nc2)
  z2[-(1:nback)]=w.forecast
  
  # fit a high order AR model for backcast
  AR.phi2=ar.yw(z2[-(1:nback)],order.max=AR.p,aic=FALSE,demean=FALSE)$ar
  for (i in nback:1) {z2[i]=sum(AR.phi2*z2[(i+1):(i+AR.p)])}

  # Add GGB factor to the z2 vector, to get GARMA realization (include forecast/update) of length narma2
  C2=gegenb.wge(u,d,nc2)
  x.forecast=c()
  for (i in 1:narma2) {x.forecast[i]=sum(C2[1:(nback+i-1)]*z2[(nback+i):2])} 
  if (forecast.step<0){MSE.AR=mean((x.ARfore[(narma+forecast.step+1):narma]-x[(narma+forecast.step+1):narma])^2)
                       MSE.GARMA=mean((x.forecast[(narma+forecast.step+1):narma]-x[(narma+forecast.step+1):narma])^2)
                      }
}else if (k==2)
{           # transform z sequence to remove the first Gegenbauer factor
            # and derive a new sequence of w, which can be fit into ARMA model
            C1=gegenb.wge(u[1],-d[1],nc)
            w1=c()
            for (i in 1:narma) {w1[i]=sum(C1[1:(nback+i-1)]*z[(nback+i):2])}
            w1=w1-mean(w1)

            # backcasting again
            w1=c(rep(0,nback),w1)
            AR.phi2=ar.yw(w1[-(1:nback)],order.max=AR.p,demean=FALSE,aic=FALSE)$ar
            for (i in nback:1) {w1[i]=sum(AR.phi2*w1[(i+1):(i+AR.p)])}
 
            C2=gegenb.wge(u[2],-d[2],nc)
            # transform w1 sequence to remove the second Gegenbauer factor
            w=c()
            for (i in 1:narma) {w[i]=sum(C2[1:(nback+i-1)]*w1[(nback+i):2])}
            w=w-mean(w)           
            w.forecast=w
  # if forecast.step>0, forecast forecast.step steps in the future
  # if forecast.step<0, update the last forecast.step steps with fitted AR model
  # meanwhile fit x with an AR model and compare the forecast made by AR and GARMA models
  x.ARest=ar.yw(x,aic=TRUE,order.max=8)
  x.ARfore=x
  if (forecast.step<0) 
  {for (i in (narma+forecast.step+1):narma) {w.forecast[i]=sum(w.forecast[(i-1):(i-p)]*phi)
                                             x.ARfore[i]=sum(x.ARfore[(i-1):(i-x.ARest$order)]*x.ARest$ar)
                                            }
  }
  if (forecast.step>0) 
  {
   w.forecast=c(w.forecast,predict(arima(w.forecast,order=c(p,0,q)),n.ahead=forecast.step)$pred)
   x.ARfore=c(x.ARfore,predict(x.ARest,n.ahead=forecast.step)$pred)
  }
  # so now the w.forecast either still has narma points (if last forecast.step steps updated)
  # or it has narma+forecast.step points (if forecast.step steps in future are forecasted)
  narma2=length(w.forecast)
  nc2=nback+narma2
  z2=rep(0,nc2)
  z2[-(1:nback)]=w.forecast
  
  # fit a high order AR model for backcast
  AR.phi2=ar.yw(z2[-(1:nback)],order.max=AR.p,aic=FALSE,demean=FALSE)$ar
  for (i in nback:1) {z2[i]=sum(AR.phi2*z2[(i+1):(i+AR.p)])}

  # Add the second GGB factor to the z2 vector, to get GARMA realization (include forecast/update) of length narma2
  C22=gegenb.wge(u[2],d[2],nc2)
  x.forecast1=c()
  for (i in 1:narma2) {x.forecast1[i]=sum(C22[1:(nback+i-1)]*z2[(nback+i):2])} 
  
  z22=rep(0,nc2)
  z22[-(1:nback)]=x.forecast1
  # fit a high order AR model
  AR.phi3=ar.yw(z22[-(1:nback)],order.max=AR.p,aic=FALSE,demean=FALSE)$ar
  for (i in nback:1) {z22[i]=sum(AR.phi3*z22[(i+1):(i+AR.p)])}
 
  C12=gegenb.wge(u[1],d[1],nc2)
  x.forecast=c()
  for (i in 1:narma2) {x.forecast[i]=sum(C12[1:(nback+i-1)]*z22[(nback+i):2])}

  if (forecast.step<0){MSE.AR=mean((x.ARfore[(narma+forecast.step+1):narma]-x[(narma+forecast.step+1):narma])^2)
                       MSE.GARMA=mean((x.forecast[(narma+forecast.step+1):narma]-x[(narma+forecast.step+1):narma])^2)
                      }
}else stop(paste("\n For now, only one Gegenbauer factor is allowed.\n"))

if (plot==TRUE)
{
if(forecast.step<0) {
t=1:length(x)
maxy=max(x)
miny=min(x)
                      plot(t,x,type="o",lty="solid",cex.lab=1,cex.main=1,cex.sub=1,pch=16,cex=.7,
                       xlab="Time",ylab="Realization/Forecast",main="Realization (Solid), GARMA Forecast (Dashed), AR Forecast (Dotted)",ylim=c(miny,maxy))
                       points(seq(narma2+forecast.step+1,narma2),mean(x)+x.forecast[(narma2+forecast.step+1):narma2],type="o",lty="dashed",pch=1,cex=.6)
                       points(seq(narma2+forecast.step+1,narma2),x.ARfore[(narma2+forecast.step+1):narma2],type="o",lty="dotted",pch=2,cex=.6)
n=length(x)
n1=n+forecast.step
n2=n1+1
dum.gf=rep(0,n)
dum.gf[n1]=x[n1]
dum.gf[n2]=mean(x)+x.forecast[n2]
dum.arf=rep(0,n)
dum.arf[n1]=x[n1]
dum.arf[n2]=x.ARfore[n2]
                       points(seq(n1,n2),dum.gf[n1:n2],type="o",lty="dashed",pch=1,cex=.7)
                       points(seq(n1,n2),dum.arf[n1:n2],type="o",lty="dotted",pch=2,cex=.7)  
                      }
 if (forecast.step>0) {
dum.gf=rep(0,n+1)
dum.gf[n]=x[n]
dum.gf[n+1]=mean(x)+x.forecast[narma+1]
dum.arf=rep(0,n+1)
dum.arf[n]=x[n]
dum.arf[n+1]=x.ARfore[narma:narma+1]
maxy=max(x)
miny=min(x)
t=1:n
plot(t,x,type="o",lty="solid",cex.lab=1,cex.main=1,cex.sub=1,pch=16,ylim=c(miny,maxy),xlim=c(1,narma2),
                       xlab="Time",ylab="Realization/Forecast",main="Realization (Solid), GARMA Forecast (Dashed), AR Forecast (Dotted)")
                       points(seq(narma+1,narma2),mean(x)+x.forecast[(narma+1):narma2],type="o",lty="dashed",pch=1,cex=.7)
                       points(seq(narma+1,narma2),x.ARfore[(narma+1):narma2],type="o",lty="dotted",pch=2,cex=.7)       
np1=n+1
t2=n:np1
                       points(seq(n,np1),dum.gf[n:np1],type="o",lty="dashed",pch=1,cex=.7)
                       points(seq(n,np1),dum.arf[n:np1],type="o",lty="dotted",pch=2,cex=.7)                                 
}
}
if(forecast.step<0){ar.fore=x.ARfore[(narma2+forecast.step+1):narma2]
garma.fore=mean(x)+x.forecast[(narma2+forecast.step+1):narma2]}
if(forecast.step>0){ar.fore=x.ARfore[(narma2-forecast.step+1):narma2]
garma.fore=mean(x)+x.forecast[(narma2-forecast.step+1):narma2]}
#cat('narma2,forecast.step',narma2,forecast.step,'\n')
#list(ARfit.order=x.ARest$order,
#     ARfore=mean(x)+x.ARfore[(narma2+forecast.step*(forecast.step<0)+1):(narma2+forecast.step*(forecast.step>0))],
#     GARMA.fore=mean(x)+x.forecast[(narma2+forecast.step*(forecast.step<0)+1):(narma2+forecast.step*(forecast.step>0))])
list(ar.fit.order=x.ARest$order,ar.fore=ar.fore,garma.fore=garma.fore)
}
