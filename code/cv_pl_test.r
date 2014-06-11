# Test particle filters on local level model
source("pf_mc_functions.r")
source("pf_functions.r")

# Set data and graphics path
dpath = "../data/"
gpath = "../graphs/"

cv_pl <- function(lambda, np, nruns=1, n.sim = 1, filt = c("pl","plrb","kd","rm"), L = 1, alpha = 0.05, progress = FALSE)
{
  rnorm(1)
  
  # Load simulated data
  load(paste(dpath,"rw_sim.rdata",sep=""))
  y = mysims[[n.sim]]$y
  nt = dim(y)[2]
  L = ceiling(L*nt) # set number of lag time points for rm pf
  
  # Run particle filters and calculate quantiles
  out = state.quant = theta.quant = list()
  for(i in 1:nruns)
  {
    out[[i]] = state.quant[[i]] = theta.quant[[i]] = list()
    if('rm' %in% filt)
    {
      length(out[[i]]) = length(state.quant[[i]]) = length(theta.quant[[i]]) = length(filt) + length(L) - 1
      names(out[[i]]) = names(state.quant[[i]]) = names(theta.quant[[i]]) = c(paste(filt[-which(filt == 'rm')]),paste('rm',L,sep=""))
    } else {
      length(out[[i]]) = length(state.quant[[i]]) = length(theta.quant[[i]]) = length(filt)
      names(out[[i]]) = names(state.quant[[i]]) = names(theta.quant[[i]]) = filt
    }
    if('pl' %in% filt)
    {
      source("pl_cv_functions.r")
      source("pl.r")
      dlpred.pl = function(y,x,suff.x,theta) dlpred.ll(y,x,suff.x,theta,lambda=lambda)
      revo.pl = function(y,x,suff.x,theta) revo.ll(y,x,suff.x,theta,lambda=lambda)
      smap.theta.pl = function(suff.theta, y, x) smap.theta.ll(suff.theta, y, x, lambda=lambda)
      smap.state.pl = function(suff.x, y, theta) smap.state.ll(suff.x, y, theta, lambda=lambda)
      out[[i]]$pl = pl(y, dlpred.pl, revo.pl, function(j) rprior.ll(lambda=lambda), rmove.ll, smap.theta.pl, smap.state.pl, np, progress=progress, method="stratified",nonuniformity="ess", threshold=0.8, log=F)
      state.quant[[i]]$pl = pf.quantile(out[[i]]$pl$state, out[[i]]$pl$weight, function(x, param=1) x, c(alpha/2,1-alpha/2))
#      theta.quant[[i]]$pl = pf.quantile(1/out[[i]]$pl$theta, out[[i]]$pl$weight, function(x, param=1) x, c(alpha/2,1-alpha/2))
      theta.quant[[i]]$pl = pf.mix.quantile(list(out[[i]]$pl$suff.theta), out[[i]]$pl$weight, list(function(x, o, w) sum(w*pgamma(x,o[1,],o[2,]))))
    }
    if('plrb' %in% filt)
    {
      source("pl_cv_functions.r")
      source("pl.r")
      dlpred.rb = function(y,x,suff.x,theta) dlpred.ll(y,x,suff.x,theta,lambda=lambda,rb=T)
      revo.rb = function(y,x,suff.x,theta) revo.ll(y,x,suff.x,theta,lambda=lambda,rb=T)
      smap.theta.rb = function(suff.theta, y, x) smap.theta.ll(suff.theta, y, x, lambda=lambda)
      smap.state.rb = function(suff.x, y, theta) smap.state.ll(suff.x, y, theta, lambda=lambda)
      out[[i]]$plrb = pl(y, dlpred.rb, revo.rb, function(j) rprior.ll(lambda=lambda), rmove.ll, smap.theta.rb, smap.state.rb, np, progress=progress, method="stratified", nonuniformity="none", log=F)
#       state.quant[[i]]$plrb = pf.quantile(out[[i]]$plrb$state, out[[i]]$plrb$weight, function(x, param=1) x, c(alpha/2,1-alpha/2))
#       theta.quant[[i]]$plrb = pf.quantile(1/out[[i]]$plrb$theta, out[[i]]$plrb$weight, function(x, param=1) x, c(alpha/2,1-alpha/2))
      state.quant[[i]]$plrb = pf.mix.quantile(list(out[[i]]$plrb$suff.state), out[[i]]$plrb$weight, list(function(x, o, w) sum(w*pnorm(x,o[1,],sqrt(o[2,])))))
      theta.quant[[i]]$plrb = pf.mix.quantile(list(out[[i]]$plrb$suff.theta), out[[i]]$plrb$weight, list(function(x, o, w) sum(w*pgamma(x,o[1,],o[2,]))))
    }
    if('kd' %in% filt)
    {
      source("kd_cv_functions.r")
      source("kd_pf.r")
      dllik.cvkd = function(y, x, theta) dllik.kd(y, x, theta, lambda=lambda)
      pstate.cvkd = function(x, theta) pstate.kd(x, theta, lambda=lambda)
      revo.cvkd = function(x, theta) revo.kd(x, theta, lambda=lambda)
      rprior.cvkd = function(j) rprior.kd()
      out[[i]]$kd = kd_pf(y, dllik.cvkd, pstate.cvkd, revo.cvkd, rprior.cvkd, np, progress=progress, method="stratified", nonuniformity="ess", threshold=0.8, log=F)
      state.quant[[i]]$kd = pf.quantile(out[[i]]$kd$state, out[[i]]$kd$weight, function(x, param=1) x, c(alpha/2,1-alpha/2))
      theta.quant[[i]]$kd = pf.quantile(1/exp(out[[i]]$kd$theta), out[[i]]$kd$weight, function(x, param=1) x, c(alpha/2,1-alpha/2))
    } 
    if('rm' %in% filt)
    {
      source("rm_cv_functions.r")
      source("rm_pf.r")
      mydlm = list(F=1,G=1,V=1,W=lambda,m0=0,C0=1)
      revo.rm = function(x, theta) revo.cv(x, theta, w = lambda)
      for(k in 1:length(L))
      {
        rmove <- function(y, x, theta) rm_mcmc(y, x, theta, a0=1, b0=1, mydlm, L[k], 1)
        ks = which(names(out[[i]]) == paste('rm',L[k],sep=""))
        out[[i]][[ks]] = rm_pf(mysims[[n.sim]]$y, dllik.cv, revo.rm, function(j) rprior.cv(), rmove, np, progress=progress, store.filter = TRUE, method="stratified", nonuniformity="ess", threshold=0.8, log=FALSE)
        state.quant[[i]][[ks]] = pf.quantile(out[[i]][[ks]]$state, out[[i]][[ks]]$weight, function(x, param=1) x, c(alpha/2,1-alpha/2))
        theta.quant[[i]][[ks]] = pf.quantile(1/out[[i]][[ks]]$theta, out[[i]][[ks]]$weight, function(x, param=1) x, c(alpha/2,1-alpha/2))
      }
    }
    print(paste(np,n.sim,i,lambda))
  }
  
  # Calculate log marginal likelihoods
  lmarglik = sapply(out, function(x) sapply(x, function(w) pf.lmarglik(w, is.null(w$increment))))
  colnames(lmarglik) = 1:nruns
  pf.out = list(lmarglik=lmarglik,state.quant=state.quant,theta.quant=theta.quant)
  save(pf.out, file=paste(dpath,"cv_pl-",lambda,"-",np,"-",nruns,"-",n.sim,".rdata",sep=""))
}

set.seed(93)
require(plyr)
require(doMC)
registerDoMC()
mydata = expand.grid(lambda = c(.5,1,2), np = c(100,500,1000,5000), nruns=20, n.sim = 1, progress=F, stringsAsFactors = FALSE)
m_ply(mydata, function(lambda,np,nruns,n.sim,L,progress) cv_pl(lambda,np,nruns,n.sim,filt = c('pl','plrb','kd','rm'),L=c(1,.25),progress=progress), .parallel = TRUE)
