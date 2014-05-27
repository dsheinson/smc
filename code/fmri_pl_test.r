# Test particle filters on local level model
source("pf_mc_functions.r")
source("pf_functions.r")
source("dlm_sim.r")
source("pl_functions.r")
source("pl.r")

# Set data and graphics path
dpath = "../data/"
gpath = "../graphs/"

fmri_pl <- function(N, mod.sim, n, n.sim, nruns, mod.est, np, alpha = 0.05, progress = FALSE)
{
  rnorm(1)
  
  # Load simulated data
  load(paste(dpath,"dlm_ar_sim-",N,"-",mod.sim,"-2.rdata",sep=""))
  mysim = get(paste(mod.sim,"_dat",sep=""))[[1]][[n]][[n.sim]]
  nt = dim(mysim$y)[2]
  
  # Create functions to run particle learning
  if(mod.est == "M101")
  {
    FF = mysim$true.params$U[1,2,]
  } else if(mod.est == "M011") {
    FF = rep(1,nt)
  } else if(mod.est == "M101s") {
    FF = 1 + mysim$true.params$U[1,2,]
  } else {stop("mod.est must be 'M101','M011',or 'M101s'")}
  dlpred = function(y, x, suff.x, theta) dlpred.ar(y, x, suff.x, theta, t(mysim$true.params$U[1,,]), FF, FALSE)
  revo = function(y, x, suff.x, theta) revo.ar(y, x, suff.x, theta, t(mysim$true.params$U[1,,]), FF, FALSE)
  smap.theta = function(suff.theta, y, x.new, x.curr) smap.theta.ar(suff.theta, y, x.new, x.curr, t(mysim$true.params$U[1,,]), FF)
  smap.state = function(suff.x, y, theta) smap.state.ar(suff.x, y, theta, t(mysim$true.params$U[1,,]), FF)
  
  # Create new sim to "train" prior
  x0 = rnorm(1, 0, sqrt(mysim$true.params$W[1,1] / (1 - mysim$true.params$G[1,1]^2)))
  sim = dlm.sim(100, array(mysim$true.params$F[,,1:100], c(1,1,100)), mysim$true.params$G, mysim$true.params$V, mysim$true.params$W, x0, mysim$true.params$beta, array(mysim$true.params$U[,,1:100],c(1,2,100)))
  
  # Run pf on training data using arbitrary prior
  d = length(mysim$true.params$beta)
  p = dim(mysim$true.params$G)[1]
  true.var = c(mysim$true.params$W[1,1], mysim$true.params$V[1,1])
  b0 = true.var^3 / 5 + true.var
  a0 = b0 / true.var + 1
  prior = list(b0 = mysim$true.params$beta, B0 = 5*diag(d), phi0 = list(mysim$true.params$G[,1]), Phi0 = list(0.5*diag(p)), as0 = a0[1], bs0 = b0[1], am0 = a0[2], bm0 = b0[2], m0 = 0)
  rprior <- function(j) rprior.ar.joint(prior)
  rmove <- function(theta, suff.theta) rmove.ar.joint(theta, suff.theta, prior)
  out <- pl(sim$y, dlpred, revo, rprior, rmove, smap.theta, smap.state, np, progress)
  
  # Create new prior from training posterior
  wt.mom <- cov.wt(t(out$theta[,,101]), out$weight[,101])
  prior = rprior.convert.joint(wt.mom)
  rprior <- function(j) rprior.train.joint(prior)
  rmove <- function(theta, suff.theta) rmove.ar.joint(theta, suff.theta, prior)
  
  # Run pfs on fMRI sim
  state.quant = theta.quant = list()
  lmarglik = times = c()
  for(i in 1:nruns)
  {
    time = system.time(out <- pl(mysim$y, dlpred, revo, rprior, rmove, smap.theta, smap.state, np, progress))
    times = c(times, time)    
    
    # Calculate quantiles and log marginal likelihood
    state.quant[[i]] = pf.quantile(out$state, out$weight, function(x, param=1) x, c(alpha/2,1-alpha/2))
    theta.quant[[i]] = pf.quantile(out$theta, out$weight, function(x, param=1) x, c(alpha/2,1-alpha/2))
    lmarglik = c(lmarglik, pf.lmarglik(out))
    
    print(paste(N,mod.sim,n,n.sim,nruns,mod.est,np,alpha,time,sep="-"))
  }
  
  pf.out = list(lmarglik=lmarglik,state.quant=state.quant,theta.quant=theta.quant,prior=prior,times)
  save(pf.out, file=paste(dpath,"fmri_pl-",paste(N,mod.sim,n,n.sim,nruns,mod.est,np,alpha,sep="-"),".rdata",sep=""))
}

require(plyr)
require(doMC)
registerDoMC()
mydata = expand.grid(N=20, mod.sim = "M011", n = 15, n.sim = 1, nruns = 20, mod.est = c("M101","M011","M101s"), np = c(100,500,1000,5000), progress = FALSE, stringsAsFactors = FALSE)
m_ply(mydata, fmri_pl, .parallel = TRUE)
