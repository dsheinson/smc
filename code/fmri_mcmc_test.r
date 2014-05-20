# Set data and graphics path
dpath = "/storage/sheinson_research/"

# Test dlm and reg MCMC functions
source("dlm_mcmc_functions.r")
source("reg_mcmc_functions.r")

fmri_dlm_mcmc_test <- function(N, n, n.sim, mod, diff, dimx, n.chain, nsims, nburn, nthin, same=FALSE, progress=TRUE, print.iter=FALSE)
{
  # Load data
  load(paste(dpath,"dlm_ar_sim-",N,"-",mod,diff,"-",dimx,".rdata",sep=""))
  mysim = get(paste(mod,"_dat",sep=""))[[1]][[n]][[n.sim]]
  y = mysim$y
  
  # Set known values and get dimensions of beta and x
  psi = list(U=mysim$true.params$U,F=mysim$true.params$F,V=mysim$true.params$V)
  d = dim(psi$U)[2]
  p = dim(psi$F)[2]
  nt = dim(psi$U)[3]

  # Set priors
  if(diff == "-diff") # multiple AR components in state vector
  {
    # assume all components are AR(1)
    phi0 = list(); length(phi0) = p
    Phi0 = list(); length(Phi0) = p
    as0 = rep(1e-6, p)
    bs0 = rep(1e-6, p)
    for(i in 1:p)
    {
      phi0[[i]] = 0
      Phi0[[i]] = matrix(1e6,1,1)
    }
  } else {
    phi0 = list(rep(0,p))
    Phi0 = list(1e6*diag(p))
    as0 = 1e-6
    bs0 = 1e-6
  }
  prior = list(m0 = rep(0,p), b0 = rep(0,d), B0 = 1e6*diag(d), am0 = 1e-6, bm0 = 1e-6, phi0 = phi0, Phi0 = Phi0, as0 = as0, bs0 = bs0)
  
  # Run MCMC
  mcmc.details = list(n.sims = nsims, n.thin = nthin, n.burn = nburn)
  out = dlm.ar.mcmc(y, psi, prior, same=same, mcmc.details=mcmc.details, progress=progress, print.iter=print.iter)
  
  # Calculate MLE # assume all error components are AR(1)
  fit = lm(y[1,] ~ t(psi$U[1,,]) - 1)
  npar = (1-same)*p + same
  phi.init = as.numeric(acf(residuals(fit),type='partial',plot=FALSE)$acf)[1]
  phi.init = rep(phi.init, npar)
  sigma2m.init = summary(fit)$sigma^2 / 2
  sigma2s.init = summary(fit)$sigma^2 / 2
  sigma2s.init = rep(sigma2s.init, npar)
  if(same) same.lab = "_same" else same.lab = ""
  build <- function(par) get(paste("build_",mod,same.lab,sep=""))(par,psi$V,psi$U,p)
  fit.mle <- dlmMLE(y, c(logit(phi.init,-1,1),log(sigma2s.init),log(sigma2m.init)),build)
  fit.smooth <- dlmSmooth(dlmFilter(y, build(fit.mle$par)))
  theta.mle <- c(fit.smooth$s[nt+1,1:d],unlogit(fit.mle$par[1:npar],-1,1),exp(fit.mle$par[-(1:npar)]))
  x.mle <- t(fit.smooth$s[,-(1:d)])
  
  # Allocate results in list and save
  out.est = list(out=out,theta.mle=theta.mle,x.mle=x.mle)
  file = paste(dpath,"fmri_dlm_mcmc_test-",paste(N,n,n.sim,mod,sep="-"),paste(diff,dimx,n.chain,nsims,nburn,nthin,same,sep="-"),".rdata",sep="")
  print(file)
  save(out.est, file = file)
}

require(dlm)
require(plyr)
require(doMC)
registerDoMC()
mydata = expand.grid(N=20,n=c(1,6,11,16),n.sim=1:8,mod=c("M101","M011"),diff="",dimx=2,n.chain=1:2,nsims=11000,nburn=1000,nthin=10,same=FALSE,progress=FALSE,print.iter=FALSE)
m_ply(mydata, fmri_dlm_mcmc_test, .parallel = TRUE)

# fmri_reg_mcmc_test <- function(N, n, n.sim, mod, dimx, n.chain, nsims, nburn, nthin, progress=TRUE, print.iter=FALSE)
# {
#   # Load data
#   load(paste(dpath,"dlm_ar_sim-",N,"-",mod,"-",dimx,".rdata",sep=""))
#   mysim = get(paste(mod,"_dat",sep=""))[[1]][[n]][[n.sim]]
#   y = mysim$y[1,]
#   
#   # Set known values and get dimensions of beta and rho
#   X = t(mysim$true.params$U[1,,])
#   d = dim(X)[2]
#   p = dim(mysim$true.params$G)[1]
#   
#   # Set priors
#   prior = list(b0 = rep(0,d), B0 = 1e6*diag(d), phi0 = rep(0,p), Phi0 = 1e6*diag(p), v0 = 1e-6, d0 = 1e-6)
#   
#   # Which parameters to sample?
#   mcmc.details = list(n.sims = nsims, n.thin = nthin, n.burn = nburn)
#   out = reg.ar.mcmc(y, X, prior, mcmc.details=mcmc.details, progress=progress, print.iter=print.iter)
#   
#   # Calculate MLEs and set initial values
#   fit = arima(y, order = c(p,0,0), xreg = X, include.mean = FALSE)
#   theta.mle = c(fit$coef[(p+1):(p+d)],fit$coef[1:p],fit$sigma2)
#   
#   # Allocate results in list and save
#   out.est = list(out=out,theta.mle=theta.mle)
#   file =  paste(dpath,"fmri_reg_mcmc_test-",paste(N,n,n.sim,mod,dimx,n.chain,nsims,nburn,nthin,sep="-"),".rdata",sep="")
#   print(file)
#   save(out.est, file = file)
# }
# 
# mydata = expand.grid(N=20,n=6,n.sim=1:20,mod=c("M010","M020"),dimx=2,n.chain=1:3,nsims=11000,nburn=1000,nthin=10,progress=FALSE,print.iter=FALSE)
# m_ply(mydata, fmri_reg_mcmc_test, .parallel = TRUE)