plotts.wge = function (x,Presentation = 0,ylab = "",xlab = "Time", title = "Realization") 
{
  cex.labs <- c(0.9, 0.8, 0.9)
  numrows <- 1
  numcols <- 1
  par(mfrow = c(numrows, numcols), mar = c(3.8, 2.5, 1, 1))
  n = length(x)
  t = 1:n
  if (Presentation == 0) {
    if (class(x) != "ts") {
      if (n <= 200) {
        plot(t, x, type = "o", cex = 0.5, pch = 16, cex.lab = 0.75, 
             cex.axis = 0.75, lwd = 0.75, xlab = xlab, 
             ylab = ylab, col = "blue", main = title)
      }
      else if (n > 200) {
        plot(t, x, type = "l", cex = 0.5, pch = 16, cex.lab = 0.75, 
             cex.axis = 0.75, lwd = 0.75, xlab = xlab, 
             ylab = ylab, col = "blue", main = title)
      }
    }
    if (class(x) == "ts") {
      if (n <= 200) {
        plot(x, type = "o", cex = 0.5, pch = 16, cex.lab = 0.75, 
             cex.axis = 0.75, lwd = 0.75, xlab = xlab, 
             ylab = ylab, col = "blue", main = title)
      }
      else if (n > 200) {
        plot(x, type = "l", cex = 0.5, pch = 16, cex.lab = 0.75, 
             cex.axis = 0.75, lwd = 0.75, xlab = xlab, 
             ylab = ylab, col = "blue", main = title)
      }
    }
  }
  else if (Presentation == 1) {
    if (class(x) != "ts") {
      df = data.frame(x = t, y = x)
      df %>% ggplot(aes(x = x, y = y)) + geom_point() + 
        geom_line() + ggtitle(title) + xlab(xlab) + 
        ylab(ylab)
    }
    else if (class(x) == "ts") {
      df = data.frame(x = zoo::as.yearmon(time(x)), y = x)
      df %>% ggplot(aes(x = x, y = y)) + geom_point() + 
        geom_line() + ggtitle(title) + xlab(xlab) + 
        ylab(ylab)
    }
  }
}

