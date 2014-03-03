source("fmri_mcmc_functions.r")

# Set data and graphics path
dpath = "../data/"
gpath = "../graphs/"

fmri_mcmc_test <- function(n.sim, mod, n.chain, x, beta, sigma2m, phi, sigma2s, progress, print.iter)
{
  # Load data
  load(paste(dpath,"dlm_ar_sim-20-",mod,".rdata",sep=""))
  y = mysims[[n.sim]]$y

  # Set known values
  V = 1
  U = mysims[[n.sim]]$true.params$U
  F = mysims[[n.sim]]$true.params$F
  psi = list(V=V,U=U,F=F)

  # Set initial values to truth and set prior on phi
  p = dim(mysims[[n.sim]]$true.params$G)[1]
  initial = list(x=mysims[[1]]$x, theta = list(beta = mysims[[1]]$true.params$beta, sigma2m = mysims[[1]]$true.params$V[1,1], phi = mysims[[1]]$true.params$G[,1], sigma2s = mysims[[1]]$true.params$W[1,1]))
  prior = list(m0 = rep(0,p), b0 = c(0,0), B0 = 1e6*diag(2), am0 = 1e-6, bm0 = 1e-6, phi0 = rep(0,p), Phi0 = 1e6*diag(p), as0 = 1e-6, bs0 = 1e-6)

  steps = c('x','beta','sigma2m','phi','sigma2s')
  params.est <- which(as.logical(c(x,beta,sigma2m,phi,sigma2s)))
  steps = steps[params.est]
  mcmc.details = list(n.sims = 110, n.thin = 1, n.burn = 10)
  out = fmri_mcmc(y, psi, prior, initial, mcmc.details, steps)
  
  save(out, file = paste(dpath,"fmri_mcmc_test-",paste(mod,n.chain,beta,sigma2m,phi,sigma2s,sep="-"),sep=""))
}