gen.sigplusnoise.wge=function (n, b0 = 0, b1 = 0, coef = c(0, 0), freq = c(0, 0), 
    psi = c(0, 0), phi = 0, vara = 1, plot = "TRUE", sn = 0) 
{
    if (sn > 0) {
        set.seed(sn)
    }
    std = sqrt(vara)
    t = 1:n
    zt = gen.arma.wge(n = n, phi = phi, theta = 0, vara = vara, 
        sn = sn)
    x = b0 + b1 * t + coef[1] * cos(2 * pi * freq[1] * t + psi[1]) + 
        coef[2] * cos(2 * pi * freq[2] * t + psi[2]) + zt
    numrows <- 1
    numcols <- 1
    fig.width <- 5.5
    fig.height <- 4.5
    cex.labs <- c(0.8, 0.7, 0.9)
    par(mfrow = c(numrows, numcols), mar = c(3.8, 2.5, 1, 1))
    if (plot == "TRUE") {
        plot(t, x, type = "l")
    }
    return(x)
}
#

