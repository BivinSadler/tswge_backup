#
wv.wge<-function(x)
{
cex.labs <- c(.7,.7,.7)
freq = seq(0, 0.5, length = 101)
input=x
xtime = (1:length(input))
input <- hilbert.wge(input)
N <- nrow <- length(freq)
ncol <- length(input)
t <- 1:ncol
tfr <- matrix(0, nrow, ncol)
for(i in 1:ncol) {
ti <- t[i]
taumax <- min(ti - 1, ncol - ti, round(N/2) - 1)
tau <- - taumax:taumax
indices <- (N + tau) %% N + 1
tfr[indices, i] <- input[ti + tau] * Conj(input[ti - tau])
tau <- round(N/2) + 1
if((ti <= ncol - tau) && (ti >= tau + 1))
tfr[tau + 1, i] <- input[ti + tau] * Conj(input[ti - tau])
}
tfr <- apply(tfr, 2, fft)
tfr <- Re(tfr)
returnval <- list(tfr = tfr, freq = freq, t = t)
numrows <- 2
numcols <- 2
fig.width <- 5
fig.height <- 2.5
cex.labs <- c(.8,.7,.8)
par(mfrow=c(numrows,numcols),mar=c(3.8,2.5,1,1))
gcol<-c(1.0,.9,.8,.7,.6,.5,.4,.3,.2,.1,0) 
image(xtime,freq,t(tfr),xaxt='n',yaxt='n',xlab='',ylab='',main='',cex.axis=1,cex.lab=1,col=gray(gcol),xlim=c(0,length(input)))
axis(side=1,cex.axis=1,mgp=c(3,0.15,0),tcl=-.3,at=NULL);
axis(side=2,las=1,cex.axis=1,mgp=c(3,.25,0),tcl=-.3)
mtext(side=1,cex=1,text='Time',line=1.3) 
mtext(side=2,cex=1,text='Frequency',line=1.3) 
mtext(side=1,cex=1,text='',line=2.3) 
#cat('tfr,',tfr,'\n')
#cat('t(tfr)',t*(tfr),'\n')

return(returnval)
}
#
#
#
