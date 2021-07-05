plotts.wge = function (x,style = 0,ylab = "",xlab = "Time", main = "Realization", col = "black", text_size= 12, lwd = .75, cex = .5, cex.lab = .75, cex.axis = .75) 
{
  #cex.labs <- c(0.9, 0.8, 0.9)
  numrows <- 1
  numcols <- 1
  #par(mfrow = c(numrows, numcols), mar = c(3.8, 2.5, 1, 1))
  n = length(x)
  t = 1:n
  if (style == 0) {
    if (class(x) != "ts") {
      if (n <= 200) {
        plot(t, x, type = "o", cex = cex, pch = 16, cex.lab = cex.lab, 
             cex.axis = cex.axis, lwd = lwd, xlab = xlab, 
             ylab = ylab, col = col, main = main, las = 1)
      }
      else if (n > 200) {
        plot(t, x, type = "l", cex = cex, pch = 16, cex.lab = cex.lab, 
             cex.axis = cex.axis, lwd = lwd, xlab = xlab, 
             ylab = ylab, col = col, main = main, las = 1)
      }
    }
    if (class(x) == "ts") {
      if (n <= 200) {
        plot(x, type = "o", cex = cex, pch = 16, cex.lab = cex.lab, 
             cex.axis = cex.axis, lwd = lwd, xlab = xlab, 
             ylab = ylab, col = col, main = main)
      }
      else if (n > 200) {
        plot(x, type = "l", cex = cex, pch = 16, cex.lab = cex.lab, 
             cex.axis = cex.axis, lwd = lwd, xlab = xlab, 
             ylab = ylab, col = col, main = main)
      }
    }
  }
  else if (style == 1) {
    if (class(x) != "ts") {
      df = data.frame(x = t, y = x)
      df %>% ggplot(aes(x = x, y = y)) + geom_point() + 
        geom_line(color = col, size = lwd) + ggtitle(main) + xlab(xlab) + 
        ylab(ylab) + theme(axis.text=element_text(size=text_size))
    }
    else if (class(x) == "ts") {
      df = data.frame(x = zoo::as.yearmon(time(x)), y = x)
      df %>% ggplot(aes(x = x, y = y)) + geom_point() + 
        geom_line(color = col, size = lwd) + ggtitle(main) + xlab(xlab) + 
        ylab(ylab) + theme(axis.text=element_text(size=text_size))
    }
  }
}

# Good reference URL for TS plotting in ggplot2: https://www.r-graph-gallery.com/279-plotting-time-series-with-ggplot2.html
