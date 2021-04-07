aic5.wge=function(x,p=0:5, q=0:2, type='aic') {
pmax=max(p)
pmin=min(p)
qmax=max(q)
qmin=min(q)
nr=(pmax-pmin+1)*(qmax-qmin+1)
aval <- matrix(0, nrow=nr,ncol=3)
mytype=type
# 
cat('---------WORKING... PLEASE WAIT...','\n')
cat('\n')
cat('\n')
#
indx=0 

for(ip in pmin:pmax)
for(iq in qmin:qmax)
{
  {
   indx<-indx+1
#cat('iteration= ',ip,iq,indx,' ','\n')
  ret<-try(aic.wge(x,p=ip,q=iq,type=mytype), silent=TRUE)
  if (is.list(ret)==TRUE){
     # put returned value into aval
     #print('returned=',ret,'\n')

     aval[indx,] <- c(ret$p,ret$q,ret$value)
    } 
    else  {
       cat('Error in aic calculation at',ip,iq,'\n')
      aval[indx,]<-c(ip,iq,999999)
    }
#    cat("debug ",ip,iq,unlist(ret), "\n") 
}
}


#
# 
# ------------------------SORT THE RESULTS
# ------------------------DISPLAY THE 5 SMALLEST AIC
# 

     #cat('Mytype=',class(mytype),'\n')
     dat <- data.frame(aval)
     sorted_aval <-dat[order(dat[,3],decreasing=F),]

     cat('Five Smallest Values of ',c(mytype),'\n')
     # 
    if(mytype=='aic') {colnames(sorted_aval)=c('   p','   q','       aic')}
    if(mytype=='aicc') {colnames(sorted_aval)=c('   p','   q','       aicc')}
    if(mytype=='bic') {colnames(sorted_aval)=c('   p','   q','       bic')}
     #aic.wgeaval
#     sorted_aval[1:5,,]
print(sorted_aval[1:5,],row.names=FALSE)
}
