kalman.miss.wge=function(y,start,gam0,F,gamV,Gtmiss,gamW){
n=length(y)
sqgamV=sqrt(gamV)
sqgamW=sqrt(gamW)
chgamW=chol(gamW)
yy=y
result=astsa::Ksmooth1(num=n, y=yy, A=Gtmiss, mu0=start, Sigma0=gam0, Phi=F, Ups=0,Gam=0,cQ=sqgamV, cR=chgamW,input=0)
pfs=cbind(yy,result$xp,result$Pp,result$xf,result$Pf,result$xs,result$Ps)
colnames(pfs)=c("Data","Prediction","Var_Predict","Filter","Var_Filter","Smooth","Var_Smooth")
return(pfs)
}
