# Test particle filters on local level model
source("pf_mc_functions.r")
source("pf_functions.r")

# Set data and graphics path
dpath = "../data/"
gpath = "../graphs/"

cv_pl <- function(lambda, np, nruns=1, n.sim = 1, filt = c("pl","pl-rb","kd","rm"), alpha = 0.05, burn = 0, progress = FALSE)
{
  rnorm(1)
  
  # Load simulated data
  load(paste(dpath,"rw_sim.rdata",sep=""))
  y = mysims[[n.sim]]$y
  
  # Run particle filters and calculate quantiles
  out = state.quant = theta.quant = list()
  for(i in 1:nruns)
  {
    out[[i]] = state.quant[[i]] = theta.quant[[i]] = list()
    if('pl' %in% filt)
    {
      source("pl_functions.r")
      source("pl.r")
      dlpred = function(y,x,suff.x,theta) dlpred.ll(y,x,suff.x,theta,lambda=lambda,rb=F)
      revo.pl = function(y,x,suff.x,theta) revo.ll(y,x,suff.x,theta,lambda=lambda,rb=F)
      smap.theta.pl = function(suff.theta, y, x.new, x.curr) smap.theta.ll(suff.theta, y, x.new, x.curr, lambda=lambda)
      smap.state.pl = function(suff.x, y, theta) smap.state.ll(suff.x, y, theta, lambda=lambda)
      out[[i]]$pl = pl(y, dlpred, revo.pl, function(j) rprior.ll(), rmove.ll, smap.theta.pl, smap.state.pl, np, progress=progress, method="stratified",nonuniformity="ess", threshold=0.8, log=F)
      state.quant[[i]]$pl = pf.quantile(out[[i]]$pl$state, out[[i]]$pl$weight, function(x, param=1) x, c(alpha/2,1-alpha/2))
      theta.quant[[i]]$pl = pf.quantile(1/out[[i]]$pl$theta, out[[i]]$pl$weight, function(x, param=1) x, c(alpha/2,1-alpha/2))
    }
    if('pl-rb' %in% filt)
    {
      source("pl_functions.r")
      source("pl.r")
      dlpred = function(y,x,suff.x,theta) dlpred.ll(y,x,suff.x,theta,lambda=lambda,rb=T)
      revo.pl = function(y,x,suff.x,theta) revo.ll(y,x,suff.x,theta,lambda=lambda,rb=T)
      smap.theta.pl = function(suff.theta, y, x.new, x.curr) smap.theta.ll(suff.theta, y, x.new, x.curr, lambda=lambda)
      smap.state.pl = function(suff.x, y, theta) smap.state.ll(suff.x, y, theta, lambda=lambda)
      out[[i]]$plrb = pl(y, dlpred, revo.pl, function(j) rprior.ll(), rmove.ll, smap.theta.pl, smap.state.pl, np, progress=progress, method="stratified",nonuniformity="ess", threshold=0.8, log=F)
      state.quant[[i]]$plrb = pf.quantile(out[[i]]$plrb$state, out[[i]]$pl$weight, function(x, param=1) x, c(alpha/2,1-alpha/2))
      theta.quant[[i]]$plrb = pf.quantile(1/out[[i]]$plrb$theta, out[[i]]$pl$weight, function(x, param=1) x, c(alpha/2,1-alpha/2))
    }
    if('kd' %in% filt)
    {
      source("kd_cv_functions.r")
      source("kd_pf.r")
      dllik.cv = function(y, x, theta) dllik.kd(y, x, theta, lambda=lambda)
      pstate.cv = function(x, theta) pstate.kd(x, theta, lambda=lambda)
      revo.cv = function(x, theta) revo.kd(x, theta, lambda=lambda)
      out[[i]]$kd = kd_pf(y, dllik.cv, pstate.cv, revo.cv, rprior.kd, np, progress=progress, method="stratified",nonuniformity="ess", threshold=0.8,log=F)
      state.quant[[i]]$kd = pf.quantile(out[[i]]$kd$state, out[[i]]$kd$weight, function(x, param=1) x, c(alpha/2,1-alpha/2))
      theta.quant[[i]]$kd = pf.quantile(1/exp(out[[i]]$kd$theta), out[[i]]$kd$weight, function(x, param=1) x, c(alpha/2,1-alpha/2))
    } 
    if('rm' %in% filt)
    {
      source("rm_cv_functions.r")
      source("rm_pf.r")
      mydlm = list(F=1,G=1,V=1,W=lambda,m0=0,C0=1)
      revo.rm = function(x, theta) revo(x, theta, w = lambda)
      rmove <- function(y, x, theta) rm_mcmc(y, x, theta, a0=1, b0=1, mydlm, 1)
      out[[i]]$rm = rm_pf(mysims[[n.sim]]$y, dllik, revo.rm, rprior, rmove, np, progress=progress, store.filter = TRUE, method="stratified", nonuniformity="ess", threshold=0.8, log=FALSE)
      state.quant[[i]]$rm = pf.quantile(out[[i]]$rm$state, out[[i]]$rm$weight, function(x, param=1) x, c(alpha/2,1-alpha/2))
      theta.quant[[i]]$rm = pf.quantile(1/out[[i]]$rm$theta, out[[i]]$rm$weight, function(x, param=1) x, c(alpha/2,1-alpha/2))
    }
    print(paste(np,n.sim,i,lambda))
  }
 
  # Calculate log marginal likelihoods
  lmarglik = sapply(out, function(x) sapply(x, function(w) pf.lmarglik(w)))
  colnames(lmarglik) = 1:nruns
  pf.out = list(lmarglik=lmarglik,state.quant=state.quant,theta.quant=theta.quant)
  save(pf.out, file=paste(dpath,"cv_pl-",lambda,"-",np,"-",nruns,"-",n.sim,".rdata",sep=""))
}

set.seed(89)
require(plyr)
require(doMC)
registerDoMC()
mydata = expand.grid(lambda = c(.5,1,2), np = c(100,500,1000,5000), nruns=20, n.sim = 1, stringsAsFactors = FALSE)
m_ply(mydata, function(lambda,np,nruns,n.sim) cv_pl(lambda,np,nruns,n.sim,filt=c('pl','pl-rb','kd','rm')), .parallel = TRUE)
