kalman.wge=function(y,start,gam0,F,gamV,G,gamW){
n=length(y)
sqgamV=sqrt(gamV)
sqgamW=sqrt(gamW)
yy=y
result=astsa::Ksmooth0(num=n, y=yy, A=G, mu0=start, Sigma0=gam0, Phi=F, cQ=sqgamV, cR=sqgamW)
pfs=cbind(yy,result$xp,result$Pp,result$xf,result$Pf,result$xs,result$Ps)
colnames(pfs)=c("Data","Prediction","Var_Predict","Filter","Var_Filter","Smooth","Var_Smooth")
return(pfs)
}

