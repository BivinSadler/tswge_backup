est.farma.wge=function(x,low.d,high.d,inc.d,p.max,nback=500)
{

# the number of Gegenbauer factors in the proposed model
k=1
low.d=low.d/2
up.d=high.d
up.d=up.d/2
# the number of observations in the input realization 
narma=length(x)
#
#
#Function ar.burg2.wge
ar.burg2.wge=function(x,order.max=5)
{
 x2=x-mean(x)
 n=length(x2)
 AIC=matrix( 0,ncol=3,nrow=(1+order.max) )

 AIC[1,1]=0
 AIC[1,2]=sum(x2^2)/n
 AIC[1,3]=log(AIC[1,2])+2*(0+1)/n
if(order.max==0){
chosen=c(AIC[1,1],AIC[1,2],AIC[1,3])}
if(order.max>0){
 for (i in 1:order.max)
 {
  burg.result=ar.burg(x,aic=FALSE,order.max=i)
  AIC[1+i,1]=i
  AIC[1+i,2]=burg.result$var.pred
  AIC[1+i,3]=log(burg.result$var.pred)+2*(i+1)/n
 }
 bestrow=which(AIC[,3]==min(AIC[,3]))
 chosen=as.vector(AIC[bestrow,])
}
 
 list(AR.order=chosen[1],sigma2=chosen[2],AIC=chosen[3])
}
#
#End ar.burg2.wge

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

# construct the grids
d1=seq(low.d[1],up.d[1],inc.d[1])
d.n=length(d1)
u1=1
u.n=1
# when there is only one Gegenbauer factor
# grid search on each u,d combination
  garmaAIC=matrix( 0, ncol=5, nrow=d.n) 
  garmaAIC[,1]=rep(u1, each=d.n)
  garmaAIC[,2]=rep(d1, u.n)

  w=matrix( 0, ncol=narma, nrow=(u.n*d.n) )

          for (i in 1:(u.n*d.n))
           {#3#
            # transform z sequence to remove Gegenbauer factors
            # and derive a new sequence of w, which can be fit into ARMA model
            C1=gegenb.wge(garmaAIC[i,1],-garmaAIC[i,2],nc)
            for (j in 1:narma) {w[i,j]=sum(C1[1:(nback+j-1)]*z[(nback+j):2])}
            w[i,]=w[i,]-mean(w[i,])
                
            # fit sequence of w with an AR(p) model
		
		garmaAIC[i,3]=ar.burg2.wge(w[i,],order.max=p.max)$AR.order
		garmaAIC[i,4]=ar.burg2.wge(w[i,],order.max=p.max)$AIC
		garmaAIC[i,5]=garmaAIC[i,4]+2*2*(garmaAIC[i,2]!=0)/narma
           }#3#
  winner=which(garmaAIC[,5]==min(garmaAIC[,5]))
  winner.final=winner[1]     # in case there are more than one winner candidate #

  d1=garmaAIC[winner.final,2]
 d=2*d1 
  p.optim=garmaAIC[winner.final,3]
  AR.AIC=garmaAIC[winner.final,4]
  GARMA.AIC=garmaAIC[winner.final,5]				 

	if (p.optim==0) { AR.coef=0
				sigma2=mean(w^2)
                      } else{
  				     model.select=ar.burg(w[winner.final,], aic=FALSE, order.max=p.optim)
				     AR.coef=model.select$ar
				     sigma2=model.select$var.pred}

 # export the result model estimates
list(d=d,phi=AR.coef,vara=sigma2,aic=GARMA.AIC)
}  
