aic.burg.wge=function (x, p = 1:5, type = "aic") 
{
k=0    
x = x - mean(x)
    aic = 99999
    bic = 99999
    aicc = 99999
    type1=type
    for (j in p)  {
        w = est.ar.wge(x,p=j,method='burg',factor=FALSE)
                   phi = w$phi
   
          res = backcast.wge(x, phi = phi,theta=0, n.back = 50)
        avar = 0
        n = length(x)
        for (i in 1:n) {
            avar = avar + res[i] * res[i]
        }
        avar = avar/n
        tempaic = log(avar) + 2 * (j + k + 1)/n
        tempbic = log(avar) + (j + k + 1) * log(n)/n
        tempaicc = log(avar) + (n + j + k + 1)/(n - j - k - 3)
        if (type1 == "aic") {
            if (tempaic < aic) {
                aic = tempaic
                j_aic = j
                 phi_aic = phi
                avar_aic = avar
            }
        }
        if (type1 == "aicc") {
            if (tempaicc < aicc) {
                aicc = tempaicc
                j_aicc = j
                phi_aicc = phi
                avar_aicc = avar
            }
        }
        if (type1== "bic") {
            if (tempbic < bic) {
                bic = tempbic
                j_bic = j
                phi_bic = phi
                avar_bic = avar
            }
        }
    }
    if (type1 == "aic") {
        out1 = list(type = type1, value = aic, p = j_aic,  
            phi = phi_aic,  vara = avar_aic)
        return(out1)
    }
    if (type1 == "aicc") {
        out1 = list(type = type1, value = aicc, p = j_aicc,  
            phi = phi_aicc, vara = avar_aicc)
        return(out1)
    }
    if (type1 == "bic") {
        out1 = list(type = type1, value = bic, p = j_bic,  
            phi = phi_bic, vara = avar_bic)
        return(out1)
    }
}
