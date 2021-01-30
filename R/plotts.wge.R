<<<<<<< HEAD
plotts.wge = function (x, yLabel = "", Presentation1 = 0, title = "Realization", xLabel = "Time") 
{
  cex.labs <- c(0.9, 0.8, 0.9)
  numrows <- 1
  numcols <- 1
  par(mfrow = c(numrows, numcols), mar = c(3.8, 2.5, 1, 1))
  n = length(x)
  t = 1:n
  if (gg == 0) {
    if (class(x) != "ts") {
      if (n <= 200) {
        plot(t, x, type = "o", cex = 0.5, pch = 16, cex.lab = 0.75, 
             cex.axis = 0.75, lwd = 0.75, xlab = xLabel, 
             ylab = yLabel, col = "blue", main = title)
      }
      else if (n > 200) {
        plot(t, x, type = "l", cex = 0.5, pch = 16, cex.lab = 0.75, 
             cex.axis = 0.75, lwd = 0.75, xlab = xLabel, 
             ylab = yLabel, col = "blue", main = title)
      }
    }
    if (class(x) == "ts") {
      if (n <= 200) {
        plot(x, type = "o", cex = 0.5, pch = 16, cex.lab = 0.75, 
             cex.axis = 0.75, lwd = 0.75, xlab = xLabel, 
             ylab = yLabel, col = "blue", main = title)
      }
      else if (n > 200) {
        plot(x, type = "l", cex = 0.5, pch = 16, cex.lab = 0.75, 
             cex.axis = 0.75, lwd = 0.75, xlab = xLabel, 
             ylab = yLabel, col = "blue", main = title)
      }
    }
  }
  else if (gg == 1) {
    if (class(x) != "ts") {
      df = data.frame(x = t, y = x)
      df %>% ggplot(aes(x = x, y = y)) + geom_point() + 
        geom_line() + ggtitle(title) + xlab(xLabel) + 
        ylab(yLabel)
    }
    else if (class(x) == "ts") {
      df = data.frame(x = zoo::as.yearmon(time(x)), y = x)
      df %>% ggplot(aes(x = x, y = y)) + geom_point() + 
        geom_line() + ggtitle(title) + xlab(xLabel) + 
        ylab(yLabel)
    }
  }
}

=======
plotts.wge <-function(x) 
{
#
#
#     x is a vector of length n containing the time series realization which is to be plotted 
#        as x(1), x(2), ..., x(n) where n=length(x)
#
#
cex.labs <- c(.9,.8,.9)
#
numrows <- 1
numcols <- 1
par(mfrow=c(numrows,numcols),mar=c(3.8,2.5,1,1))
#
n=length(x)
t=1:n
#
#Plot Data
#
#
if(n <=200) {plot(t,x,type='o',xaxt='n',yaxt='n',cex=0.5,pch=16,cex.lab=.75,cex.axis=.75,lwd=.75,xlab='',ylab='', col = "blue")}
if(n > 200) {plot(t,x,type='l',xaxt='n',yaxt='n',cex=0.5,pch=16,cex.lab=.75,cex.axis=.75,lwd=.75,xlab='',ylab='')}
axis(side=1,cex.axis=.9,mgp=c(3,0.15,0),tcl=-.3);
axis(side=2,las=1,cex.axis=.9,mgp=c(3,.4,0),tcl=-.3)
mtext(side=c(1,2,1),cex=cex.labs,text=c('Time','','Realization'),line=c(1,1.1,2.1))
#                 }
#     
}
>>>>>>> 2ed5d2e40d1c46a040ee8d7b26d8f2241b0d8f66
