# Test particle filters on local level model
source("pf_mc_functions.r")
source("pf_functions.r")
source("dlm_sim.r")
source("pl_fmri_functions.r")
source("pl.r")

# Set data and graphics path
dpath = "../data/"
gpath = "../graphs/"

fmri_pl <- function(N, mod.sim, n, n.sim, nruns, mod.est, np, prior, alpha = 0.05, progress = FALSE)
{
  rnorm(1)
  
  # Load simulated data
  load(paste(dpath,"dlm_ar_sim-",N,"-",mod.sim,"-2.rdata",sep=""))
  mysim = get(paste(mod.sim,"_dat",sep=""))[[1]][[n]][[n.sim]]
  nt = dim(mysim$y)[2]
  p.sim = dim(mysim$x)[1]
  d = dim(mysim$true.params$U)[2]
  
  # Create functions to run particle learning
  if(mod.est == "M101")
  {
    FF = mysim$true.params$U[1,d,]
  } else if(mod.est == "M011") {
    FF = rep(1,nt)
  } else if(mod.est == "M101s") {
    FF = 1 + mysim$true.params$U[1,d,]
  } else if(mod.est == "M021") {
    FF = rep(1,nt)
  } else {stop("mod.est must be 'M101','M011','M101s', or 'M021'")}
  dlpred = function(y, x, suff.x, theta) dlpred.ar(y, x, suff.x, theta, t(mysim$true.params$U[1,,]), FF)
  revo = function(y, x, suff.x, theta) revo.ar(y, x, suff.x, theta, t(mysim$true.params$U[1,,]), FF)
  smap.theta = function(suff.theta, y, x.new, x.curr) smap.theta.ar(suff.theta, y, x.new, x.curr, t(mysim$true.params$U[1,,]), FF)
  smap.state = function(suff.x, y, theta) smap.state.ar(suff.x, y, theta, t(mysim$true.params$U[1,,]), FF)
  
  if(missing(prior))
  {
    # Calculate MLE, standard errors, and use to create new prior
    require(dlm)
    p = sum(as.numeric(strsplit(mod.est,"")[[1]][2:3]))
    dyn = strsplit(mod.est,"")[[1]][2] == 1
    fit = dlmMLE(mysim$y, c(rep(.5,p),1,1), function(par) build.ar1(par, t(mysim$true.params$U[1,,]), dyn), hessian=T)
    s=dlmSmooth(dlmFilter(mysim$y, build.ar1(fit$par, t(mysim$true.params$U[1,,]))))
    mle = c(s$s[nt+1,1:d],fit$par[1:p],exp(fit$par[(p+1):(p+2)]))
#     cov = diag(c(fit$par[1:p],exp(fit$par[-(1:p)]))) %*% solve(fit$hessian) %*% diag(c(fit$par[1:p],exp(fit$par[-(1:p)]))) 
#     cov.s = dlmSvd2var(s$U.S,s$D.S)
#     prior = rprior.convert(list(cov=bdiag(list(cov.s[[nt+1]][(1:d),(1:d)],cov)),center=mle), d)
#     prior$m0 = rep(0,p); prior$C0 = cov.s[[nt+1]][d+1,d+1]*diag(p)
    prior = list(b0 = mle[1:d], B0 = 10*diag)
    rprior <- function(j) rprior.ar(prior)
    rmove <- function(suff.theta) rmove.ar(suff.theta, prior)
  } else {
    rprior <- function(j) rprior.ar(prior)
    rmove <- function(suff.theta) rmove.ar(suff.theta, prior)
  }
  
  # Run pfs on fMRI sim
  state.quant = theta.quant = list()
  lmarglik = times = c()
  for(i in 1:nruns)
  {
    time = system.time(out <- pl(mysim$y, dlpred, revo, rprior, rmove, smap.theta, smap.state, np, progress, method="stratified", nonuniformity="ess", threshold=0.8, log=F))
    times = c(times, time[3])
    
    # Calculate quantiles and log marginal likelihood
    state.quant[[i]] = pf.quantile(out$state, out$weight, function(x, param=1) x, c(alpha/2,.5,1-alpha/2))
    theta.quant[[i]] = pf.quantile(out$theta, out$weight, function(x, param=1) x, c(alpha/2,.5,1-alpha/2))
    lmarglik = c(lmarglik, pf.lmarglik(out))
    
    print(paste(mod.sim,n,n.sim,i,mod.est,np,alpha,round(time[3],2),sep="-"))
  }
  
  pf.out = list(lmarglik=lmarglik,state.quant=state.quant,theta.quant=theta.quant,prior=prior,times=times)
  save(pf.out, file=paste(dpath,"fmri_pl-",paste(mod.sim,n,n.sim,nruns,mod.est,np,alpha,sep="-"),".rdata",sep=""))
}

require(plyr)
require(doMC)
registerDoMC()
mydata = expand.grid(N=20, mod.sim = c("M101","M011","M021"), n = 19, n.sim = 3, nruns = 5, mod.est = c("M101","M011","M021"), np = c(100,500,1000), progress = F, stringsAsFactors = FALSE)
m_ply(mydata, function(N,mod.sim,n,n.sim,nruns,mod.est,np,progress) fmri_pl(N, mod.sim, n, n.sim, nruns, mod.est, np, progress=progress), .parallel=T)
