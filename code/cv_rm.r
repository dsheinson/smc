source("rm_pf.r")
source("rm_cv_functions.r")
source("pf_mc_functions.r")
source("pf_functions.r")

# Set data path
dpath = "../data/"

# Load simulated data
load(paste(dpath,"rw_sim.rdata",sep=""))

# Function to run pfs
pf_lmarglik <- function(n.sim, n.run, W, np, alpha = 0.05, progress = FALSE)
{
  rnorm(1)
  revo1 <- function(x, theta) revo(x, theta, w = W)
  mydlm = list(F = 1, G = 1, V = 1, W = W, m0 = 0, C0 = 1)
  rmove <- function(y, x, theta) rm_mcmc(y, x, theta, a0=1, b0=1, mydlm, 1)
  out = rm_pf(mysims[[n.sim]]$y, dllik, revo1, rprior, rmove, np, progress=progress, method="stratified", nonuniformity="ess", threshold=0.8, log=FALSE)
  
  # Calculate log-marginal likelihood
  lmarglik = pf.lmarglik(out)
  
  # Calculate credible intervals of states and precision
  state.quant = pf.quantile(out$state, out$weight, function(x, param=1) x, c(alpha/2,1-alpha/2))
  theta.quant = pf.quantile(1/out$theta, out$weight, function(x, param=1) x, c(alpha/2,1-alpha/2))
  
  pf.out = list(lmarglik=lmarglik, state.quant=state.quant, theta.quant=theta.quant)
  cat(n.sim, n.run, W, np,'\n')
  save(pf.out, file=paste(dpath,"cv_pf-",paste(n.sim, n.run, W, np, alpha,sep="-"),".rdata",sep=""))
}

# Apply pf_lmarglik for many pf runs with different number of particles and values of W
data1 = expand.grid(n.run=1:20, W=c(0.1,0.5,1,2,3), np=c(500,1000,5000,10000), n.sim=1, stringsAsFactors=FALSE)
data2 = expand.grid(n.run=1, W=c(0.1,0.5,1,2,3), np=c(500,1000,5000,10000), n.sim=2:20, stringsAsFactors=FALSE)
mydata = rbind(data1,data2)
require(doMC)
registerDoMC()
require(plyr)
m_ply(.data = mydata, .fun = pf_lmarglik, .parallel = TRUE)
