source("rm_pf.r")
source("pf_functions.r")
source("pf_mc_functions.r")
source("rm_fmri_functions.r")

# Set data and graphics path
dpath = "../data/"
gpath = "../graphs/"

# Set seed
set.seed(54)

# Function to run pfs
fmri_rm <- function(N, mod.sim, dimx, n, n.sim, n.run, mod.est, np, alpha = 0.05, progress = TRUE)
{
  # Load simulated data
  load(paste(dpath,"dlm_ar_sim-",N,"-",mod.sim,"-",dimx,".rdata",sep=""))
  mysim = get(paste(mod.sim,"_dat",sep=""))[[1]][[n]][[n.sim]]
  
  rnorm(1)
  dllik <- function(y, x, theta) get(paste("dllik.",mod.est,sep=""))(y, x, theta, mysim$true.params$U[1,2,])
  rmove <- function(y, x, theta) rmcmc(y, x, theta, mysim$true.params$U[1,2,], mod.est)
  out = rm_pf(mysim$y, dllik, revo, rprior, rmove, np, progress=progress, method="stratified", nonuniformity="ess", threshold=0.8, log=FALSE)
  
  # Calculate log-marginal likelihood
  lmarglik = pf.lmarglik(out)
  
  # Calculate credible intervals of states and precision
  state.quant = pf.quantile(out$state, out$weight, function(x, param=1) x, c(alpha/2,1-alpha/2))
  theta.quant = pf.quantile(out$theta, out$weight, function(x, param=1) x, c(alpha/2,1-alpha/2))
  
  pf.out = list(lmarglik=lmarglik, state.quant=state.quant, theta.quant=theta.quant)
  file = paste(dpath,"fmri_rm-",paste(N,mod.sim,dimx,n,n.sim,n.run,mod.est,np,sep="-"),".rdata",sep="")
  print(file)
  save(pf.out, file=file)
}

require(plyr)
require(doMC)
registerDoMC()
mydata = expand.grid(N = 20, mod.est = c("M101","M011"), dimx = 2, n = 6, n.sim = 1:20, n.run = 1:4, mod.sim = c("M101","M011"), np = 5000, progress = FALSE)
m_ply(.dat = mydata, .fun = fmri_rm, .parallel = TRUE)