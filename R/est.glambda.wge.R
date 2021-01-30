# est.glambda.wge lets the user specify ranges of candidate lambda and shift values,
# it builds a matrix, each combination of lambda and shift values is in one row.
# according to this pair of lambda and shift values, the Q stat is calculated,
# choose the lambda and shift with the smallest Q-value.
#--------------------------------------------------------------------------------

est.glambda.wge<-function(data,lambda.range=c(0,1),offset.range=c(0,100))
{
#------------------------------------------------------------
# function picknk.wge
#------------------------------------------------------------
picknk.wge=function(data,offset,lambda)
{
	a=trans.to.dual.wge(data,offset=offset,lambda=lambda,plot=FALSE)
	return(getnk.acf.wge(a$intY)$mse)    #a$intY is the values of unequally spaced X(t_k)'s,#
                                       #i.e. the equally spaced dual Y_k's#
}
#--------------------------------------------------------------
# end function picknk.wge
#---------------------------------------------------------------
#
#------------------------------------------------------------
# function getnk.acf.wge
#------------------------------------------------------------
getnk.acf.wge=function(data,lag0=0)  #split a data set into two halves, then calculate their own acf#
                                 #and calculate the mse of the two sets of acf's#
{
	n = length(data)
	if(lag0<1) lag0 = floor(n/8)
	acf1 = acf(data[1:floor(n/2)], lag.max = lag0, plot = F)[[1]]
	acf2 = acf(data[(floor(n/2) + 1):n], lag.max = lag0, plot = F)[[1]]
      acf1.5 = acf(data[(floor(n/4) + 1):floor(3*n/4)], lag.max = lag0, plot = F, type="covariance")[[1]]
      #mse = sum((acf1 - acf1.5)^2)+sum((acf2-acf1.5)^2) 
      mse = sum( (acf1 - acf2)^2 )
	list(mse = mse, acf1 = acf1, acf2 = acf2)
}
#------------------------------------------------------------
# end function getnk.acf.wge
#------------------------------------------------------------
#
#------------------------------------------------------------
# function getper.p.wge
#------------------------------------------------------------
getper.p.wge=function(data,lag=3,method="max")
{
	n=length(data)
	m=floor(n/lag)
	a=NULL
	i=1
	n1=0
	
	while(n1<n) {                                      
		       data1=data[(n1+1):min(n,(n1+lag+2))]
		       if(method=="max") m2=which(data1==max(data1))[1]
		       if(method=="min") m2=which(data1==min(data1))[1]

		       if((m2>1)&(m2<lag+2)&((m2+n1)!=n)) a=c(a,(m2+n1))
		       i=i+1
		       n1=(i-1)*lag
	            }
	return(a)
}
#------------------------------------------------------------
# end function getper.p.wge
#------------------------------------------------------------
#
#
	Interp="Linear"
                lambdaRange=lambda.range
                shiftRange=offset.range
                if(mode(data)=="list") data<-as.numeric(unlist(data))
	n<-length(data)
	ip<-getper.p.wge(data,lag=3,method="max")
	lambda<-seq(lambdaRange[1],lambdaRange[2],by=0.1)
	shift<-seq(shiftRange[1],shiftRange[2])

	Qslope<-matrix(0,length(lambda),3)
	Qslope[,1] <- lambda
	shiftslope <- matrix(0,length(shift),3)

	for(j in 1:length(lambda)) {
		                      shiftslope[,1] <- lambda[j]
		                      shiftslope[,2] <- shift
			                shiftslope[,3] <- sapply(shift, picknk.wge, data = data, lambda = lambda[j])
			                # If there is not a unique offset that attains the minimum Q-value, 
			                # choose the first offset in the range (i.e., this will usually be the smallest).
			                bestrow <- which(shiftslope[,3]==min(shiftslope[,3]))
			                if(length(bestrow)>1) bestrow <- bestrow[1]
			                Qslope[j,] <- shiftslope[bestrow,]
		                     }

	dimnames(Qslope)[[2]]<-c("lambda","offset","Q")
	Qslope<-data.frame(Qslope)
	bestQ<-Qslope[Qslope$Q==min(Qslope$Q),]

	list(Q=Qslope,bestlambda=bestQ[,1],bestshift=bestQ[,2])
}
