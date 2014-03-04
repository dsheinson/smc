source("dlm_ar_functions.r")
source("rm_pf.r")
source("rm_cv_functions.r")
source("pf_mc_functions.r")
source("pf_functions.r")

# Set data path
dpath = "../data/"

# Load simulated data
load(paste(dpath,"rw_sim.rdata",sep=""))

# Function to run pfs
pf_lmarglik <- function(n.sim, np, W, label, alpha = 0.05, progress = FALSE)
{
  rnorm(1)
  a0 = b0 = 1
  mydlm = list(F = 1, G = 1, V = 1, W = W, m0 = 0, C0 = 1)
  dllik1 = function(y, x, theta) dllik(y, x, sigma2 = theta)
  revo1 = function(x, theta) revo(x, sigma2 = theta, W = W)
  rprior1 <- function() rprior(a0,b0)
  rmove <- function(y, x, theta) rm_mcmc(y, x, theta, a0, b0, mydlm, 1)
  out = rm_pf(mysims[[n.sim]]$y, dllik1, revo1, rprior1, rmove, np, progress=progress, method="stratified", nonuniformity="ess", threshold=0.8, log=FALSE)
  
  # Calculate log-marginal likelihood
  lmarglik = pf.lmarglik(out)
  
  # Calculate credible intervals of states and precision
  state.quant = pf.quantile(out$state, out$weight, function(x, param=1) x, c(alpha/2,1-alpha/2))
  theta.quant = pf.quantile(1/out$theta, out$weight, function(x, param=1) x, c(alpha/2,1-alpha/2))
  
  pf.out = list(lmarglik=lmarglik, state.quant=state.quant, theta.quant=theta.quant)
  cat(n.sim,np,W,label,'\n')
  save(pf.out, file=paste(dpath,"cv_pf-",n.sim,"-",np,"-",W,"-",label,"-",alpha,".rdata",sep=""))
}

# Apply pf_lmarglik for many pf runs with different number of particles and values of W
data1 = expand.grid(n.sim=rep(c(1,2,3),c(20,20,20)),np=10000, W=2, stringsAsFactors=FALSE)
data1 = data.frame(data1, label = 1:dim(data1)[1])
data2 = expand.grid(n.sim=4:20,np=10000, W=c(0.5,1,2), stringsAsFactors=FALSE)
data2 = data.frame(data2, label = (dim(data1)[1] + 1:dim(data2)[1]))
mydata = rbind(data1,data2)
require(doMC)
registerDoMC()
require(plyr)
set.seed(61)
m_ply(.data = mydata, .fun = pf_lmarglik, .parallel = TRUE)
