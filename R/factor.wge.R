#
#  Factor Table Routine
# 
#
#  Print Factor Table in file test1.txt in Documents
#
#
#
factor.wge=function(phi)
{
#
# phi is a vector containing the coefficients (AR or MA) defining the polynomial to be factored
#
mphi=format(round(phi, 4), trim=T, nsmall = 4)
cat('\n')
cat('Coefficients of Original polynomial: ','\n')
cat(mphi,'\n')
cat('\n')
cat('Factor                 Roots                Abs Recip    System Freq','\n')
mphi=-phi
one=c(1)
coef=c(one,mphi)   # these are the actual coefficients of the pth order polynomial
root=polyroot(coef)    # with signs in characteristic equation
	root1 <- sort(root)
	temp1 <- rep(0, 5)
	# fac should have one row for each linear or quadratic factor
	nfactors <- 0.5*(length(root1)+sum(abs(Im(root1)) <= 10^(-5)))
	fac <- matrix(0, nfactors, 5)
	rootindex <- 1
	for(i in 1:nfactors) {
		if(round(root1[rootindex+1], 5) == round(Conj(
			root1[rootindex]), 5) && abs(Im(root1[rootindex])) > 10^(-5)) {
			temp1[4] <- 1/Mod(root1[rootindex])
			temp1[5] <- abs(Arg(root1[rootindex])/(2 * pi))
			temp1[1] <- -2 * Re(1/root1[rootindex])
			temp1[2] <- 1/Re(root1[rootindex] * root1[rootindex+1])
                        temp1[3] <- root1[rootindex]
			fac[i,] <- temp1
			rootindex <- rootindex+2
		}
		else {
			temp1[4] <- Re(1/Mod(root1[rootindex]))
			temp1[5] <- abs(Arg(root1[rootindex])/(2 * pi))
			temp1[1] <-  - Re(1/root1[rootindex])
			temp1[2] <- 0
                        temp1[3] <- root1[rootindex]
			fac[i,] <- temp1
			rootindex <- rootindex+1
		}

	}
	fac<-data.frame(fac)
dimnames(fac)[[2]] <- c("Fac1", "Fac2", "Root", "Abs_Recip", "Freq")
fac<-fac[order(-fac$Abs_Recip),]
#
for(i in 1:nfactors) {
fc1<--as.double(fac[i,1])
fc2<--as.double(fac[i,2])
rootr<-Re(fac[i,3])
rooti<-abs(Im(fac[i,3]))
absrec<-as.double(fac[i,4])
freq<-as.double(fac[i,5])
afc14=format(round(abs(fc1), 4), trim=T, nsmall=4)
fc24=format(round(fc2, 4), trim=T, nsmall = 4)
mfc24=format(round(-fc2, 4), trim=T, nsmall = 4)
rootr4=format(round(rootr, 4), trim=T, nsmall = 4)
mrootr4=format(round(-rootr, 4), trim=T, nsmall = 4)
rooti4=format(round(rooti, 4), trim=T, nsmall = 4)
absrec4=format(round(absrec, 4), trim=T, nsmall = 4)
freq4=format(round(freq, 4), trim=T, nsmall = 4)
format(round(freq, 4), nsmall = 4,trim=T)
if(fc2 != 0) {
if(fc1 < 0 & fc2 < 0 & rootr4 > 0)   cat(sep="","1+", afc14,"B+",mfc24,"B^2","    ",rootr4,"+-",rooti4,"i","      ",absrec4,"       ",freq4,'\n')
if(fc1 < 0 & fc2 < 0 & rootr4 < 0)   cat(sep="","1+", afc14,"B+",mfc24,"B^2","   ",rootr4,"+-",rooti4,"i","      ",absrec4,"       ",freq4,'\n')
if(fc1 >= 0 & fc2 < 0)  cat(sep="","1-", afc14,"B+",mfc24,"B^2","    ",rootr4,"+-",rooti4,"i","      ",absrec4,"       ",freq4,'\n')
             }
             else { 
if(fc1<0) cat(sep="","1+",afc14,"B","             ",rootr4,"               ",absrec4,"       ",freq4,'\n')
if(fc1>0) cat(sep="","1-",afc14,"B","              ",rootr4,"               ",absrec4,"       ",freq4,'\n')  
              }
                  }
cat(' ','\n')
cat(' ','\n')
}








