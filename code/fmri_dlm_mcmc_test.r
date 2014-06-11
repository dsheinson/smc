# Set data and graphics path
dpath = "../data/"

# Test dlm and reg MCMC functions
source("dlm_mcmc_functions.r")
source("pl_fmri_functions.r")

fmri_dlm_mcmc_test <- function(N, n, n.sim, mod.sim, mod.est, n.chain, nsims, nburn, nthin, same=FALSE, prior.type="vague", progress=TRUE, print.iter=FALSE)
{
  # Load data
  load(paste(dpath,"dlm_ar_sim-",N,"-",mod.sim,"-2.rdata",sep=""))
  mysim = get(paste(mod.sim,"_dat",sep=""))[[1]][[n]][[n.sim]]
  y = mysim$y
  
  # Set known values and get dimensions of beta and x
  nt = dim(mysim$true.params$U)[3]
  d = dim(mysim$true.params$U)[2]
  p.sim = dim(mysim$true.params$F)[2]
  p = sum(as.numeric(strsplit(mod.est,"")[[1]][2:3]))
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
  psi = list(U=mysim$true.params$U,F=array(rbind(FF,rep(0,p-1)),c(1,p,nt)),V=matrix(1))

  # Set priors
  if(prior.type == "vague")
  {
    phi0 = list(rep(0,p))
    Phi0 = list(1e6*diag(p))
    as0 = 1e-6
    bs0 = 1e-6
    prior = list(m0 = rep(0,p), b0 = rep(0,d), B0 = 1e6*diag(d), am0 = 1e-6, bm0 = 1e-6, phi0 = phi0, Phi0 = Phi0, as0 = as0, bs0 = bs0)
  } else if(prior.type == "rm"){
    load(paste(dpath,"fmri_rm-",paste(N,mod,dimx,n,n.sim,1,mod,100,1,5,sep="-"),".rdata",sep=""))
    prior = pf.out$prior
  } else if(prior.type == "mle") {
    require(dlm)
    dyn = strsplit(mod.est,"")[[1]][2] == 1
    fit = dlmMLE(y, c(rep(.5,p),1,1), function(par) build.ar1(par, t(psi$U[1,,]), dyn), hessian=T)
    s=dlmSmooth(dlmFilter(y, build.ar1(fit$par, t(psi$U[1,,]))))
    mle = c(s$s[nt+1,1:d],fit$par[1:p],exp(fit$par[(p+1):(p+2)]))
    cov = diag(c(fit$par[1:p],exp(fit$par[-(1:p)]))) %*% solve(fit$hessian) %*% diag(c(fit$par[1:p],exp(fit$par[-(1:p)]))) 
    cov.s = dlmSvd2var(s$U.S,s$D.S)
    prior = rprior.convert(list(cov=bdiag(list(cov.s[[nt+1]][(1:d),(1:d)],cov)),center=mle), d)
    prior$m0 = rep(0,p); prior$C0 = cov.s[[nt+1]][d+1,d+1]*diag(p)
  }
  

  # Run MCMC
  mcmc.details = list(n.sims = nsims, n.thin = nthin, n.burn = nburn)
  out = dlm.ar.mcmc(y, psi, prior, same=same, mcmc.details=mcmc.details, progress=progress, print.iter=print.iter)
  
  # Calculate MLE # assume all error components are AR(1)
  dyn = strsplit(mod.est,"")[[1]][2] == 1
  fit = dlmMLE(y, c(rep(.5,p),1,1), function(par) build.ar1(par, t(psi$U[1,,]), dyn), hessian=T)
  s=dlmSmooth(dlmFilter(y, build.ar1(fit$par, t(psi$U[1,,]))))
  theta.mle = c(s$s[nt+1,1:d],fit$par[1:p],exp(fit$par[(p+1):(p+2)]))
  x.mle = t(s$s[,-(1:d)])
  
  # Allocate results in list and save
  out.est = list(out=out,theta.mle=theta.mle,x.mle=x.mle,prior=prior)
  file = paste(dpath,"fmri_dlm_mcmc_test-",paste(N,n,n.sim,mod.sim,mod.est,prior.type,n.chain,nsims,nburn,nthin,same,sep="-"),".rdata",sep="")
  print(file)
  save(out.est, file = file)
}

require(dlm)
require(plyr)
mydata = expand.grid(N=20,n=19,n.sim=3,mod.sim=c("M101","M011","M021"),mod.est=c("M101","M011","M021"),n.chain=1,nsims=50,nburn=8,nthin=4,same=FALSE,prior.type='mle',progress=FALSE,print.iter=FALSE,stringsAsFactors=F)

require(doMC)
registerDoMC()
#mydata = expand.grid(N=20,n=19,n.sim=3,mod.sim=c("M101","M011","M021"),mod.est=c("M101","M011","M021"),n.chain=1:3,nsims=1100,nburn=100,nthin=1,same=FALSE,prior.type='mle',progress=FALSE,print.iter=FALSE)
m_ply(mydata, fmri_dlm_mcmc_test, .parallel = TRUE)