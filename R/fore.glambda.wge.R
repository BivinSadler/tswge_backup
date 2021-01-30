fore.glambda.wge=function (data.orig, lambda = 0, offset = 60, phi = 0, h = 0, n.ahead = 10, lastn = TRUE, plot = TRUE) 
{
#--------------------------------------------------------------
#  reinterpol.wge
#--------------------------------------------------------------
reinterpol.wge<-function(data,start,h,lambda,int.start,length,IntMethod="Linear")
{
  if(lambda == 0) {
    m <- log(start)/log(h)
    ht <- h^seq(m, (m + length(data)-1))
  }
  else {
    m <- (start^lambda - 1)/(h * lambda)
    kk <- seq(m, (m + length(data)-1))
    ht <- (kk * h * lambda + 1)^(1/lambda)
  }
  tt <- seq(int.start, (int.start + length-1))
  if(IntMethod=="Linear") intData <- approx(ht, data, xout = tt, rule = 3)$y
  return(intData)
}
#--------------------------------------------------------------
#  end reinterpol.wge
#--------------------------------------------------------------

    if (lambda == 0) {
        lambda = 0.001
    }
    shift = offset
    IntMethod = "Linear"
    theta = 0
    forecast.step = n.ahead
    lags = n.ahead
   nahead=n.ahead    
   n <- length(data.orig)
    dualData = trans.to.dual.wge(x = data.orig, lambda = lambda, 
        offset = offset, h = 0, plot = FALSE)
    dualData.int <- dualData$intY
#
#   
 if (lastn == "FALSE") {
        if (lambda != 0) {
            len1 <- ceiling(((n + shift + lags)^lambda - (n + 
                shift)^lambda)/(dualData$h * lambda))
        }
        else {
            len1 <- lags
        }
        int.start <- n + shift + 1
        fore <- 1
        glam.dual.fore = fore.arma.wge(dualData.int, phi = phi, 
            theta = 0, n.ahead = len1, lastn = FALSE, limits = FALSE, 
            plot = FALSE)
     
   p1 = length(phi) + 1
        start = dualData$intX[p1]
        h = dualData$h
        glam.data.fore <- reinterpol.wge(data = c(dualData.int[1:length(dualData.int)], 
            glam.dual.fore$f), start = (shift + 1), h = dualData$h, 
            int.start = int.start, lambda = lambda, length = lags, 
            IntMethod = IntMethod)
        plot.start = n - 2 * nahead + 1
        nlast = n + nahead
        data.plot = rep(0, n + nahead)
        data.plot = c(data.orig[plot.start:n], glam.data.fore)
        n.plot = 3 * nahead
        t = plot.start:nlast
#
        if(plot==TRUE){
            cex.labs=c(.9,.9,1)
            plot(t, data.plot, type = "l", xaxt = "n", yaxt = "n", 
            cex = 0.6, pch = 16, cex.lab = 0.75, cex.axis = 0.75, 
            lwd = 0.75, xlab = "", ylab = "", lty = 5, col = 4)
        axis(side = 1, cex.axis = 0.8, mgp = c(3, 0.15, 0), tcl = -0.3)
        axis(side = 2, las = 1, cex.axis = 0.8, mgp = c(3, 0.4, 
            0), tcl = -0.3)
        mtext(side = c(1, 2, 1), cex = cex.labs, text = c("Time", 
            "", "glambda forecasts: blue dashed lines;  AR forecasts: red dotted lines"), 
            line = c(0.8, 1.1, 1.8))
        n1 = n + 1
        nah = n + n.ahead
        to = plot.start:n
        points(to, data.orig[plot.start:n], type = "l")
        tff = n:nah
        ar.fit = aic.wge(data.orig, p = 0:8, q = 0:0, type = "aic")
        ar.fore = fore.arma.wge(data.orig, phi = ar.fit$phi, 
            theta = 0, n.ahead = nahead, lastn = FALSE, limits = FALSE, 
            plot = FALSE)
        ar.plot = c(data.orig[n], ar.fore$f)
        points(tff, ar.plot, "l", lty = 3, cex = 0.6, lwd = 0.4, 
            col = 2)
    }
    }
    forecast.step = n.ahead
#
#
 if (lastn == TRUE) {
        forecast.step = -forecast.step
        n = length(data.orig)
        lenl = length(dualData$intX[dualData$intX > (n + offset + 
            forecast.step)])
        int.start = n + forecast.step + shift + 1
        nahead = n.ahead
        ar.fit = aic.wge(data.orig, p = 0:4, q = 0:0, type = "aic")
        ar.fore = fore.arma.wge(data.orig, phi = ar.fit$phi, 
            theta = 0, n.ahead = nahead, lastn = TRUE, limits = FALSE, 
            plot = FALSE)
        glam.dual.fore = fore.arma.wge(dualData.int, phi = phi, 
            theta = 0, n.ahead = lenl, lastn = TRUE, limits = FALSE, 
            plot = FALSE)
        glam.data.fore = reinterpol.wge(data = c(dualData$intY[1:(length(dualData$intY) - 
            lenl)], glam.dual.fore$f), start = shift + 1, h = dualData$h, 
            int.start = int.start, lambda = lambda, length = abs(forecast.step), 
            IntMethod = IntMethod)
        nst = max(n + 5 * forecast.step, 1)
#
        if(plot==TRUE){
            cex.labs=c(.9,.9,1)
            plot(seq(nst, n), data.orig[nst:n], type = "l",  xaxt = "n", yaxt = "n", 
            cex = 0.6, pch = 16, cex.lab = 0.75, cex.axis = 0.75, 
            lwd = 0.75, xlab = "", ylab = "")       
        axis(side = 1, cex.axis = 0.8, mgp = c(3, 0.15, 0), tcl = -0.3)
        axis(side = 2, las = 1, cex.axis = 0.8, mgp = c(3, 0.4, 
            0), tcl = -0.3)
        mtext(side = c(1, 2, 1), cex = cex.labs, text = c("Time", 
            "", "glambda forecasts: blue dashed lines;  AR forecasts: red dotted lines"), 
            line = c(0.8, 1.1, 1.8))

 points(seq((n + forecast.step), n), c(data.orig[n + forecast.step], 
            ar.fore$f), type = "l", lty = 3, lwd = 0.4, col = 2)
 points(seq((n + forecast.step), n), c(data.orig[n + forecast.step+1], 
            glam.data.fore), type = "l", lty = 5, lwd = 0.4, 
            col = 4)
     }
    }
    out1 = list(f.ar = ar.fore$f, f.glam = glam.data.fore)
   return(out1)
}
