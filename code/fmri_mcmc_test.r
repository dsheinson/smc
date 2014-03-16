# Set data and graphics path
dpath = "../data/"
gpath = "../graphs/"

# Test dlm MCMC functions
source("dlm_mcmc_functions.r")

fmri_dlm_mcmc_test <- function(N, n, n.sim, mod, n.chain, nsims, nburn, nthin, x=1, beta=1, sigma2m=1, phi=1, sigma2s=1, progress=TRUE, print.iter=FALSE)
{
  require(dlm)
  
  # Load data
  load(paste(dpath,"dlm_ar_sim-",N,"-",mod,".rdata",sep=""))
  mysim = get(paste(mod,"_dat",sep=""))[[1]][[n]][[n.sim]]
  y = mysim$y

  # Set known values and get dimensions of beta and phi
  U = mysim$true.params$U
  F = mysim$true.params$F
  V = mysim$true.params$V
  psi = list(U=U,F=F,V=V)
  d = dim(U)[2]
  p = dim(F)[2]
  nt = dim(U)[3]

  # Calculate MLEs and set initial values
  fit = lm(y[1,] ~ t(U[1,,]) - 1)
  phi.init = as.numeric(acf(residuals(fit),type='partial',plot=FALSE)$acf)[1:p]
  sigma2m.init = summary(fit)$sigma^2 / 2
  sigma2s.init = summary(fit)$sigma^2 / 2
  build <- function(par) get(paste("build_",mod,sep=""))(par,V,U,U[1,2,])
  fit.mle <- dlmMLE(y, c(logit(phi.init,-1,1),log(sigma2s.init),log(sigma2m.init)),build)
  fit.smooth <- dlmSmooth(dlmFilter(y, build(fit.mle$par)))
  if(x)
  {
    x.init = t(fit.smooth$s[,(d+1):(d+p)])
    if(n.chain > 1) for(j in 1:p) x.init[j,] = rnorm(nt+1,x.init[j,],sd(x.init[j,]))
  } else {x.init = mysim$x}
  if(beta)
  {
    beta.init = fit.smooth$s[nt+1,1:d]
    if(n.chain > 1) beta.init <- rnorm(d, beta.init, sqrt(beta.init))
  } else {beta.init = mysim$true.params$beta}
  if(phi)
  {
    phi.init = unlogit(fit.mle$par[1],-1,1)
    if(n.chain > 1) phi.init = rnorm(p, phi.init, .25)
    while(!is.stationary(phi.init)) phi.init = phi.init = rnorm(p, phi.init, .25)
  } else {phi.init = mysim$true.params$G[,1]}
  if(sigma2m)
  { 
    sigma2m.init = exp(fit.mle$par[3])
    if(n.chain > 1) sigma2m.init = rnorm(1,sigma2m.init,1)
    while(sigma2m.init <= 0) sigma2m.init = rnorm(1,sigma2m.init,1)
  } else { sigma2m.init = mysim$true.params$V[1,1]}
  if(sigma2s)
  { 
    sigma2s.init = exp(fit.mle$par[2])
    if(n.chain > 1) sigma2s.init = rnorm(1,sigma2s.init,1)
    while(sigma2s.init <= 0) sigma2s.init = rnorm(1,sigma2s.init,1)
  } else { sigma2s.init = mysim$true.params$W[1,1]}
  initial = list(x=x.init, theta = list(beta = beta.init, sigma2m = sigma2m.init, phi = phi.init, sigma2s = sigma2s.init))
  
  # Set priors
  prior = list(m0 = rep(0,p), b0 = rep(0,d), B0 = 1e6*diag(d), am0 = 1e-6, bm0 = 1e-6, phi0 = rep(0,p), Phi0 = 1e6*diag(p), as0 = 1e-6, bs0 = 1e-6)

  # Which parameters to sample?
  steps = c('x','beta','sigma2m','phi','sigma2s')
  params.est <- which(as.logical(c(x,beta,sigma2m,phi,sigma2s)))
  steps = steps[params.est]
  mcmc.details = list(n.sims = nsims, n.thin = nthin, n.burn = nburn)
  out = dlm.ar.mcmc(y, psi, prior, initial, mcmc.details, steps, progress, print.iter)
  
  cat(N,n,n.sim,mod,n.chain,nsims,nburn,nthin,x,beta,sigma2m,phi,sigma2s,"\n",sep=" ")
  save(out, file = paste(dpath,"fmri_dlm_mcmc_test-",paste(N,n,n.sim,mod,n.chain,nsims,nburn,nthin,x,beta,sigma2m,phi,sigma2s,sep="-"),".rdata",sep=""))

  return(out)
}

require(plyr)
require(doMC)
registerDoMC()
mydata = expand.grid(N=20,n=6,n.sim=1:20,mod=c("M011","M101"),n.chain=1:3,nsims=11000,nburn=1000,nthin=10,progress=FALSE,print.iter=FALSE)
out.all = mlply(mydata, fmri_dlm_mcmc_test, .parallel = TRUE)

# Test reg MCMC functions
source("reg_mcmc_functions.r")

fmri_reg_mcmc_test <- function(N, n, n.sim, mod, n.chain, nsims, nburn, nthin, beta=1, phi=1, sigma2s=1, progress=TRUE, print.iter=FALSE)
{
  # Load data
  load(paste(dpath,"dlm_ar_sim-",N,"-",mod,".rdata",sep=""))
  mysim = get(paste(mod,"_dat",sep=""))[[1]][[n]][[n.sim]]
  y = mysim$y[1,]
  
  # Set known values and get dimensions of beta and rho
  X = t(mysim$true.params$U[1,,])
  d = dim(X)[2]
  p = dim(mysim$true.params$G)[1]
  
  # Calculate MLEs and set initial values
  fit = arima(y, order = c(p,0,0), xreg = X, include.mean = FALSE)
  if(beta)
  {
    beta.init <-fit$coef[(p+1):(p+d)]
    if(n.chain > 1) beta.init = rnorm(d,beta.init,sqrt(beta.init))
  } else beta.init = mysim$true.params$beta
  if(phi)
  {
    phi.init <- fit$coef[1:p]
    if(n.chain > 1) phi.init = rnorm(p,phi.init,.1)
    while(!is.stationary(phi.init)) phi.init = rnorm(p,phi.init,.1)
  } else {phi.init = mysim$true.params$G[,1]}
  if(sigma2s)
  {
    sigma2s.init <- fit$sigma2
    if(n.chain > 1) sigma2s.init = rnorm(1,sigma2s.init,1)
    while(sigma2s.init <= 0) sigma2s.init = rnorm(1,sigma2s.init,1)
  } else {sigma2s.init = mysim$true.params$W[1,1]}
  initial = list(beta = beta.init, phi = phi.init, sigma2s = sigma2s.init)
  
  # Set priors
  prior = list(b0 = rep(0,d), B0 = 1e6*diag(d), phi0 = rep(0,p), Phi0 = 1e6*diag(p), v0 = 1e-6, d0 = 1e-6)
  
  # Which parameters to sample?
  steps = c('beta','phi','sigma2s')
  params.est <- which(as.logical(c(beta,phi,sigma2s)))
  steps = steps[params.est]
  mcmc.details = list(n.sims = nsims, n.thin = nthin, n.burn = nburn)
  out = reg.ar.mcmc(y, X, prior, initial, mcmc.details, steps, progress, print.iter)
  
  cat(N,n,n.sim,mod,n.chain,nsims,nburn,nthin,beta,phi,sigma2s,"\n",sep=" ")
  save(out, file = paste(dpath,"fmri_reg_mcmc_test-",paste(N,n,n.sim,mod,n.chain,nsims,nburn,nthin,beta,phi,sigma2s,sep="-"),".rdata",sep=""))
  
  return(out)
}

require(plyr)
require(doMC)
registerDoMC()
mydata = expand.grid(N=20,n=6,n.sim=1:20,mod=c("M010","M020"),n.chain=1:3,nsims=11000,nburn=1000,nthin=10,progress=FALSE,print.iter=FALSE)
out.ar.all = mlply(mydata, fmri_reg_mcmc_test, .parallel = TRUE)
