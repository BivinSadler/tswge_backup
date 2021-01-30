backcast.wge=function(x,phi=0,theta=0,n.back=50){
# x is realization
# phi contains AR parameter (estimates)
# theta contains MA parameter(estimates)
# n.back specifies how far back we backcast (i.e. to x(-n.back))
# Output: vector of  n residuial estimates based on backcast procedure
#
#
#
n=length(x)
p=length(phi)
if(sum(phi^2)==0) {p=0}
q=length(theta)
if(sum(theta^2)==0) {q=0}
deltar=rep(0,n)
n.full=n+n.back+1
xhat=rep(0,n.full)
xbar=mean(x)
const=1
if (p > 0) {for(jp in 1:p) {const=const-phi[jp]}}
#
#
# Calculate Backcast residuals (deltar(t))
#
##
maconst=const*xbar
mpq=max(p,q)
np=n-mpq
for (ii in 1:np) {i=np-ii+1

      deltar[i]=x[i]
#cat('i,ii,deltar[i],x[i]',i,ii,deltar[i],x[i],'\n')
    if (p>0) {for (jp in 1:p) {deltar[i]=deltar[i]-phi[jp]*x[i+jp]
#              cat('i,jp,deltar',i,jp,deltar[i],'\n')
}}
#
#
#
#cat('in backcast,p,phi,q,theta',p,phi,q,theta,'\n')
    if (q>0) {for (jq in 1:q) {deltar[i]=deltar[i]+theta[jq]*deltar[i+jq]
#               cat('i,jq,deltar',i,jq,deltar[i],'\n')
}}
          deltar[i]=deltar[i]-maconst
#
#cat('i,deltar',deltar[i],'\n')
}

#
#
#

#   Backcast X(0), X(-1), ...., X(-50)
#   Call these backcasts XX(51),XX(50),......XX(1)
#       where n.back=50, npn.back=51
#       we use X(1) ... X(p) which are indexed and called xhat(52)...,xhat(52+p-1)
#       or x(npn.back+1), ..., x(npn.back+p-1)
#
#       
npn.back=n.back+1
mm=npn.back
for (i in 1:n) {
     iib=npn.back+i
     xhat[iib]=x[i]}
for (h in 1:npn.back) {
if (p > 0) {for (jp in 1:p) {xhat[mm-h+1]=xhat[mm-h+1]+phi[jp]*xhat[mm-h+jp+1]}}
if (h<=q) {for(jq in h:q) {xhat[mm-h+1]=xhat[mm-h+1]-theta[h]*deltar[jq]}}
                    xhat[mm-h+1]=xhat[mm-h+1]+maconst}
#
#
#

# Calculating Residuals after backcasting
#
#
p1=p+1
q1=q+1
p1q1=max(p1,q1)
resid=rep(0,n.full)
for (i in p1q1:n) {xhat[npn.back+i]=x[i]}

for (i in p1q1:n.full) {resid[i]=xhat[i]
   if ( p > 0) {for (jp in 1:p) {resid[i]=resid[i]-phi[jp]*xhat[i-jp]}}
   if (q > 0) {for (jq in 1:q) {resid[i]=resid[i]+theta[jq]*resid[i-jq]}}
                   resid[i]=resid[i]-maconst}

residb=rep(0,n)
for(i in 1:n) {residb[i]=resid[i+n.back+1]}
return(residb)
}
