plotts.wge=function(x, style = 0, xlab = "", ylab = "Time", 
    main = "", col = "black", text_size = 12, lwd = 0.75, 
    cex = 0.5, cex.lab = 0.75, cex.axis = 0.75, xlim = NULL, 
    ylim = NULL) 
{
    numrows <- 1
    numcols <- 1
    n = length(x)
    t = 1:n
    if (style == 0) {
        if (class(x) != "ts") {
            if (n <= 200) {
                plot(t, x, type = "o", cex = cex, pch = 16, 
                  cex.lab = cex.lab, cex.axis = cex.axis, lwd = lwd, 
                  xlab = xlab, ylab = ylab, col = col, main = main, 
                  las = 1, xlim = xlim, ylim = ylim, las = 1)
            }
            else if (n > 200) {
                plot(t, x, type = "l", cex = cex, pch = 16, 
                  cex.lab = cex.lab, cex.axis = cex.axis, lwd = lwd, 
                  xlab = xlab, ylab = ylab, col = col, main = main, 
                  las = 1, xlim = xlim, ylim = ylim, las = 1)
            }
        }
        if (class(x) == "ts") {
            if (n <= 200) {
                plot(t,x, type = "o", cex = cex, pch = 16, 
                  cex.lab = cex.lab, cex.axis = cex.axis, lwd = lwd, 
                  xlab = xlab, ylab = ylab, col = col, main = main, 
                  xlim = xlim, ylim = ylim, las = 1)
            }
            else if (n > 200) {
                plot(t,x, type = "l", cex = cex, pch = 16, 
                  cex.lab = cex.lab, cex.axis = cex.axis, lwd = lwd, 
                  xlab = xlab, ylab = ylab, col = col, main = main, 
                  xlim = xlim, ylim = ylim, las = 1)
            }
        }
    }
    else if (style == 1) {
        if (class(x) != "ts") {
            df = data.frame(x = t, y = x)
            df %>% ggplot(aes(x = x, y = y)) + geom_point() + 
                geom_line(color = col, size = lwd) + ggtitle(main) + 
                xlab(xlab) + ylab(ylab) + theme(text = element_text(size = rel(text_size))) + 
                coord_cartesian(xlim = xlim, ylim = ylim)
        }
        else if (class(x) == "ts") {
            df = data.frame(x = zoo::as.yearmon(time(x)), y = x)
            df %>% ggplot(aes(x = x, y = y)) + geom_point() + 
                geom_line(color = col, size = lwd) + ggtitle(main) + 
                xlab(xlab) + ylab(ylab) + theme(text = element_text(size = rel(text_size))) + 
                coord_cartesian(xlim = xlim, ylim = ylim)
        }
    }
}
