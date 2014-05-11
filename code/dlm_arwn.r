source("dlm_ar_functions.r")
source("dlm_sim.r")
require(dlm)

# Set data and graphics path
dpath = "../data/"
gpath = "../graphs/"

arwn.mle <- function(nt, nsims, x, beta0 = 900, beta1 = 5, phi = 0.8, sigma2s = 1, sigma2m = 1)
{
  # Create x variable
  if(x == "conv")
  {
    load(paste(dpath,"fmri-design-arwn.rdata",sep=""))
    ind = which(sapply(fmri.design, function(a) length(a$t)) == nt)
    X = fmri.design[[ind]]$X
  } else if(x == "norm") {
    X = cbind(1,rnorm(nt))
  } else if(x == "none") {
    X = cbind(1,rep(0,nt))
  }
  
  # Create storage for MLEs
  theta.mle = matrix(NA, nr = nsims, nc = 5)
  colnames(theta.mle) = c("beta0", "beta1", "phi", "sigma2s", "sigma2m")
  
  # Start loop over simulations
  for(i in 1:nsims)
  {
    # Simulate time series
    F = array(1, c(1,1,nt))
    G = matrix(phi)
    V = matrix(sigma2m)
    W = matrix(sigma2s)
    x0 = rnorm(1, 0, sqrt(sigma2s / (1 - phi^2)))
    beta = c(beta0,beta1)
    U = array(t(X), c(1, 2, nt))
    mysim = dlm.sim(nt, F, G, V, W, x0, beta, U)
    y = mysim$y
    
    # Calculate MLEs of fit to M011
    if(x == "none") U = array(U[1,1,], c(1,1,nt))
    if(x != "none") fit = lm(y[1,] ~ t(U[1,,]) - 1) else fit = lm(y[1,] ~ 1)
    phi.init = as.numeric(acf(residuals(fit),type='partial',plot=FALSE)$acf)[1]
    sigma2m.init = summary(fit)$sigma^2 / 2
    sigma2s.init = summary(fit)$sigma^2 / 2
    fit.mle <- dlmMLE(y, c(logit(phi.init,-1,1),log(sigma2s.init),log(sigma2m.init)), function(par) build_M011(par, V, U))
    fit.smooth <- dlmSmooth(dlmFilter(y, build_M011(fit.mle$par, V, U)))
    theta.mle[i,] <- c(fit.smooth$s[nt+1,1:2],unlogit(fit.mle$par[1],-1,1),exp(fit.mle$par[-1]))
 
    print(paste(nt, nsims, x, i))
  }
  
  return(theta.mle)
}

require(plyr)
require(doMC)
registerDoMC()
mydata = expand.grid(nt = c(500,1000,5000,10000), nsims = 1000, x = c("conv","norm","none"), stringsAsFactors=FALSE)
mle.all = maply(mydata, arwn.mle, .parallel = TRUE)
save(mle.all, file = paste(dpath,"dlm_arwn.rdata",sep=""))