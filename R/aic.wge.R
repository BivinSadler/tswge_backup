aic.wge=function (x, p = 0:5, q = 0:2, type = "aic") 
{
   xbar=mean(x)
   x = x - mean(x)
    aic = 99999
    bic = 99999
    aicc = 99999
    for (j in p) for (k in q) {
        w = getOption("warn")
        options(warn = -1)
        b = try(arima(x, c(j, 0, k)), silent=TRUE)
     
#       if (is.list(b)==FALSE){
  #       cat('Detect Error iteration= ',j,k,' ',is.list(b),'\n')
  #        out1 = list(type = type, value = bic, p = j_bic, q = k_bic, 
  #              phi = phi_bic, theta = theta_bic, avar = avar_bic)          
  #       return(out1)
 #       }
  
        options(warn = w)
if(is.list(b)==TRUE){
        c = as.vector(coef(b))
        if (j == 0) {
            phi = 0
        }
        else {
            phi = c[1:j]
        }
        if (k == 0) {
            theta = 0
        }
        else {
            theta = -c[(j + 1):(j + k)]
        }
        res = backcast.wge(x, phi = phi, theta = theta, n.back = 50)
        avar = 0
        n = length(x)
        for (i in 1:n) {
            avar = avar + res[i] * res[i]
        }
        avar = avar/n
        tempaic = log(avar) + 2 * (j + k + 1)/n
        tempbic = log(avar) + (j + k + 1) * log(n)/n
        tempaicc = log(avar) + (n + j + k + 1)/(n - j - k - 3)
        if (type == "aic") {
            if (tempaic < aic) {
                aic = tempaic
                j_aic = j
                k_aic = k
                phi_aic = phi
                theta_aic = theta
                avar_aic = avar
            }
        }
        if (type == "aicc") {
            if (tempaicc < aicc) {
                aicc = tempaicc
                j_aicc = j
                k_aicc = k
                phi_aicc = phi
                theta_aicc = theta
                avar_aicc = avar
            }
        }
        if (type == "bic") {
            if (tempbic < bic) {
                bic = tempbic
                j_bic = j
                k_bic = k
                phi_bic = phi
                theta_bic = theta
                avar_bic = avar
            }
        }
    }
}
    if (type == "aic") {
        out1 = list(type = type, value = aic, p = j_aic, q = k_aic, 
            phi = phi_aic, theta = theta_aic, xbar=xbar,vara = avar_aic)
        return(out1)
    }
    if (type == "aicc") {
        out1 = list(type = type, value = aicc, p = j_aicc, q = k_aicc, 
            phi = phi_aicc, theta = theta_aicc, xbar=xbar,vara = avar_aicc)
        return(out1)
    }
    if (type == "bic") {
        out1 = list(type = type, value = bic, p = j_bic, q = k_bic, 
            phi = phi_bic, theta = theta_bic, xbar=xbar,vara = avar_bic)
        return(out1)
    }
}

