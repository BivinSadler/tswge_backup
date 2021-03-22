aic5.ar.wge=function(x,p=0:5, type='aic',method='mle') {
pmax=max(p)
pmin=min(p)
#if(pmin==0) pmin=1
nr=(pmax-pmin+1)
aval <- matrix(0, nrow=nr,ncol=2)
mytype=type
mymethod=method
# 
cat('---------WORKING... PLEASE WAIT...','\n')
cat('\n')
cat('\n')
#
indx=0 

for(ip in pmin:pmax)
  {
   indx<-indx+1
#cat('iteration= ',ip,iq,indx,' ','\n')
  ret<-try(aic.ar.wge(x,p=ip,type=mytype,method=mymethod), silent=TRUE)
  if (is.list(ret)==TRUE){
     # put returned value into aval
     #print('returned=',ret,'\n')

     aval[indx,] <- c(ret$p,ret$value)
    } 
    else  {
       cat('Error in aic calculation at',ip,'\n')
      aval[indx,]<-c(ip,999999)
    }
#    cat("debug ",ipunlist(ret), "\n") 
}
#
# 
# ------------------------SORT THE RESULTS
# ------------------------DISPLAY THE 5 SMALLEST AIC
# 

     #cat('Mytype=',class(mytype),'\n')
     dat <- data.frame(aval)
     sorted_aval <-dat[order(dat[,2],decreasing=F),]

     cat('Five Smallest Values of ',c(mytype),'\n')
     cat('Method=',mymethod,'\n')
     # 
    if(mytype=='aic') {colnames(sorted_aval)=c('   p','       aic')}
    if(mytype=='aicc') {colnames(sorted_aval)=c('   p','       aicc')}
    if(mytype=='bic') {colnames(sorted_aval)=c('   p','       bic')}
     #aic.wgeaval
#     sorted_aval[1:5,,]
print(sorted_aval[1:5,],row.names=FALSE)
}
