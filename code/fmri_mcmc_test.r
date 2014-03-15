# Set data and graphics path
dpath = "../data/"
gpath = "../graphs/"

# Test dlm MCMC functions
source("dlm_mcmc_functions.r")

fmri_dlm_mcmc_test <- function(N, n, n.sim, mod, n.chain, nsims, nburn, nthin, x=1, beta=1, sigma2m=1, phi=1, sigma2s=1, progress=TRUE, print.iter=FALSE)
{
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

  # Set initial values
  if(x)
  {
    tt = dim(mysim$x)[2]
    x.init = matrix(NA,nr=p,nc=tt)
    for(j in 1:p) x.init[j,] = rnorm(tt,mysim$x[j,],sd(mysim$x[j,]))
  } else {x.init = mysim$x}
  if(beta) beta.init = rnorm(d,mysim$true.params$beta,sqrt(mysim$true.params$beta)) else beta.init = mysim$true.params$beta
  if(phi)
  {
    phi.init = rnorm(p,mysim$true.params$G[,1],sqrt(mysim$true.params$G[,1]))
    while(!is.stationary(phi.init)) phi.init = rnorm(p,mysim$true.params$G[,1],sqrt(mysim$true.params$G[,1]))
  } else {phi.init = mysim$true.params$G[,1]}
  if(sigma2m)
  { 
    sigma2m.init = rnorm(1,mysim$true.params$V[1,1],sqrt(mysim$true.params$V[1,1]))
    while(sigma2m.init <= 0) sigma2m.init = rnorm(1,mysim$true.params$V[1,1],sqrt(mysim$true.params$V[1,1]))
  } else { sigma2m.init = mysim$true.params$V[1,1]}
  if(sigma2s)
  {
    sigma2s.init = rnorm(1,mysim$true.params$W[1,1],sqrt(mysim$true.params$W[1,1]))
    while(sigma2s.init <= 0) sigma2s.init = rnorm(1,mysim$true.params$W[1,1],sqrt(mysim$true.params$W[1,1]))
  } else {sigma2s.init = mysim$true.params$W[1,1]}
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
mydata = expand.grid(N=1,n=6,n.sim=1,mod="M101",n.chain=1,nsims=100,nburn=10,nthin=9)
out.all = mlply(mydata, fmri_dlm_mcmc_test, .parallel = TRUE)
