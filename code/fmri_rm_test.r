source("rm_pf.r")
source("pf_functions.r")
source("pf_mc_functions.r")
source("rm_fmri_functions.r")
source("dlm_sim.r")

# Set data and graphics path
dpath = "../data/"
gpath = "../graphs/"

# Function to run pfs
fmri_rm <- function(N, mod.sim, dimx, n, n.sim, n.run, mod.est, np, sd.fac = 1, alpha = 0.05, progress = TRUE)
{
   # Load simulated data
  load(paste(dpath,"dlm_ar_sim-",N,"-",mod.sim,"-",dimx,".rdata",sep=""))
  mysim = get(paste(mod.sim,"_dat",sep=""))[[1]][[n]][[n.sim]]
  nt = dim(mysim$y)[2]
  
  # Create new sim to "train" prior
  x0 = rnorm(1, 0, sqrt(mysim$true.params$W[1,1] / (1 - mysim$true.params$G[1,1]^2)))
  sim = dlm.sim(100, array(mysim$true.params$F[,,1:100], c(1,1,100)), mysim$true.params$G, mysim$true.params$V, mysim$true.params$W, x0, mysim$true.params$beta, array(mysim$true.params$U[,,1:100],c(1,2,100)))
  
  # Run rm pf on training data
  rnorm(1)
  prior = list(b0 = c(900,5), B0 = 5*diag(2), phi0 = list(0.8), Phi0 = list(matrix(0.2)), as0 = 1, bs0 = 1, am0 = 1, bm0 = 1, m0 = 0)
  rprior <- function() rprior.train(prior)
  dllik <- function(y, x, theta) get(paste("dllik.",mod.est,sep=""))(y, x, theta, mysim$true.params$U[1,2,])
  rmove <- function(y, x, theta) rmcmc(y, x, theta, mysim$true.params$U[1,2,], mod.est, prior)
  out = rm_pf(sim$y, dllik, revo, rprior, rmove, np, progress=progress, method="stratified", nonuniformity="ess", threshold=0.8, log=FALSE)
  
  # Create new prior from training posterior
  wt.mom <- cov.wt(t(out$theta[,,101]), out$weight[,101])
  prior = rprior.convert(wt.mom$center, diag(wt.mom$cov), sd.fac)
  rprior <- function() rprior.train(prior)
  dllik <- function(y, x, theta) get(paste("dllik.",mod.est,sep=""))(y, x, theta, mysim$true.params$U[1,2,])
  rmove <- function(y, x, theta) rmcmc(y, x, theta, mysim$true.params$U[1,2,], mod.est, prior)
  out = rm_pf(mysim$y, dllik, revo, rprior, rmove, np, progress=progress, method="stratified", nonuniformity="ess", threshold=0.8, log=FALSE)
  
  # Calculate log-marginal likelihood
  lmarglik = pf.lmarglik(out)
  
  # Calculate credible intervals of states and precision
  state.quant = pf.quantile(out$state, out$weight, function(x, param=1) x, c(alpha/2,1-alpha/2))
  theta.quant = pf.quantile(out$theta, out$weight, function(x, param=1) x, c(alpha/2,1-alpha/2))
  
  pf.out = list(lmarglik=lmarglik, state.quant=state.quant, theta.quant=theta.quant)
  file = paste(dpath,"fmri_rm-",paste(N,mod.sim,dimx,n,n.sim,n.run,mod.est,np,sd.fac,100*alpha,sep="-"),".rdata",sep="")
  print(file)
  save(pf.out, file=file)
}

require(plyr)
require(doMC)
registerDoMC()
data0 = expand.grid(N = 20, mod.est = "M011", dimx = 2, n.sim = 8, n.run = 3, mod.sim = "M011", sd.fac = 1, n = 1, np = 100, progress = FALSE)
data1 = expand.grid(N = 20, mod.est = c("M101","M011"), dimx = 2, n.sim = 1:8, n.run = 1:3, mod.sim = c("M101","M011"), sd.fac = 1, n = c(11,16), np = 100, progress = FALSE)
data2 = expand.grid(N = 20, mod.est = c("M101","M011"), dimx = 2, n.sim = 1:8, n.run = 1:3, mod.sim = c("M101","M011"), sd.fac = c(5,7), n = 6, np = 100, progress = FALSE)
mydata = rbind(data0,data1,data2)
m_ply(.dat = mydata, .fun = fmri_rm, .parallel = TRUE)