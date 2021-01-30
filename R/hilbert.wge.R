#################################################################
### This function calculates the hilbert transformation of a given real valued signal(even length). ###
hilbert.wge<- function( input)
{
n <- length( input )
if( (n%%2)==1 )
{
input <- input[-1]
n <- length( input )
}
x1 <- fft( input ) 
x2 <- rep(0,n)
x2[1] <- x2[n/2+1] <- 1
x2[2:(n/2)] <- 2
x3 <- x1*x2
ans <- fft( x3, inverse=T )/n
return( ans )
}
######