# Set data and graphics path
dpath = "/storage/sheinson_research/"

# Test reg MCMC functions
source("reg_mcmc_functions.r")

fmri_reg_mcmc_test <- function(N, n, n.sim, mod, dimx, n.chain, nsims, nburn, nthin, progress=TRUE, print.iter=FALSE)
{
  # Load data
  load(paste(dpath,"dlm_ar_sim-",N,"-",mod,"-",dimx,".rdata",sep=""))
  mysim = get(paste(mod,"_dat",sep=""))[[1]][[n]][[n.sim]]
  y = mysim$y[1,]
  
  # Set known values and get dimensions of beta and rho
  X = t(mysim$true.params$U[1,,])
  d = dim(X)[2]
  p = dim(mysim$true.params$G)[1]
  
  # Set priors
  prior = list(b0 = rep(0,d), B0 = 1e6*diag(d), phi0 = rep(0,p), Phi0 = 1e6*diag(p), v0 = 1e-6, d0 = 1e-6)
  
  # Which parameters to sample?
  mcmc.details = list(n.sims = nsims, n.thin = nthin, n.burn = nburn)
  out = reg.ar.mcmc(y, X, prior, mcmc.details=mcmc.details, progress=progress, print.iter=print.iter)
  
  # Calculate MLEs and set initial values
  fit = arima(y, order = c(p,0,0), xreg = X, include.mean = FALSE)
  theta.mle = c(fit$coef[(p+1):(p+d)],fit$coef[1:p],fit$sigma2)
  
  # Allocate results in list and save
  out.est = list(out=out,theta.mle=theta.mle)
  file =  paste(dpath,"fmri_reg_mcmc_test-",paste(N,n,n.sim,mod,dimx,n.chain,nsims,nburn,nthin,sep="-"),".rdata",sep="")
  print(file)
  save(out.est, file = file)
}

require(plyr)
require(doMC)
registerDoMC()
mydata = expand.grid(N=20,n=6,n.sim=1:10,mod=c("M010","M020"),dimx=3,n.chain=1:3,nsims=11000,nburn=1000,nthin=1,progress=FALSE,print.iter=FALSE)
m_ply(mydata, fmri_reg_mcmc_test, .parallel = TRUE)