# Function that will draw the Unit Circle, plot a root, and calculate if it is inside or outside
# the Unit Circle

# Simply pass it the real and imaginary parts and it will plot the root and calculate 
# distance from the center

#Example: UnitCircle(.3,-.7)


unit.circle.wge = function(real = 0, imaginary = 0)
{

#Coordiate Plane and Circle

#im = mtext(expression(paste(italic("im"))), side = 2, line = -19.9, at = 1.3)
#re = mtext(expression(paste(italic("re"))), side = 3, line = -14.5, at = 1.4)
  
plot.new()
plot(-1.5:1.5, -1.5:1.5,type="n",asp = 1, xlab=expression(paste(italic("re"))),ylab=expression(paste(italic("im"))),main="Roots and the Unit Circle", cex.axis = 1.5, cex.lab = 1.5, las = 1)

draw.circle(0,0,1,border="blue",lty=2,lwd=2)
abline(h = 0)
abline(v = 0)

if(length(real) == 1){
# plot the point designated by (real,imaginary)
points(real,imaginary, pch = 19, col = "red")
lines(c(0,real), c(0,imaginary), lwd = 3, col = "red")
}

if(length(real) == 2){
  # plot the point designated by (real,imaginary)
  points(real[1],imaginary[1], pch = 19, col = "red")
  lines(c(0,real[1]), c(0,imaginary[1]), lwd = 3, col = "red")
  points(real[2],imaginary[2], pch = 19, col = "red")
  lines(c(0,real[2]), c(0,imaginary[2]), lwd = 3, col = "red")
  }


#Calc Length
length = sqrt(real^2 + imaginary^2)

#Display distance from origin and Abs Reciprocal Labels
#mtext("Length", side = 3, at = c(-1.8), line = -1)
#mtext("Abs Reciprocal", side = 3, at = c(1.5), line = -1)

#Display distance from origin and Abs Reciprocal Values
#mtext(as.character(round(length,4)), side = 3, at = c(-1.8), line = -1.9)
#mtext(as.character(round(1/abs(length),4)), side = 3, at = c(1.5), line = -1.9)

}
