source("mc_functions.r")
source("rm_pf.r")
source("rm_test_functions.r")

# Load simulated data
load("../data/dlm_sim.rdata")

# Function to run resample move particle filter and approximate log marginal likelihood of the data
pf_lmarglik <- function(nsim, np, W, corr = FALSE, progress = FALSE)
{
  a0 = b0 = 1
  rprior1 <- function() rprior(a0,b0)
  revo1 <- function(x, theta) revo(x, theta, w = W)
  rmove <- function(y, x, theta) rm_mcmc(y, x, theta, a0, b0, 1)
  if(corr) set.seed(W)
  out = rm_pf(sims$y[1,], dllik, revo1, rprior1, rmove, np, progress = progress, method="stratified", nonuniformity="ess", threshold=0.8, log=FALSE)
  lmarglik = pf.lmarglik(out)
  cat(np,nsim,W,'\n')
  write(paste(lmarglik,np,nsim,W,sep=","),"../data/mc_1sim_test.txt",append=TRUE)
  return(lmarglik)
}

# Apply pf_lmarglik for many pf runs with different number of particles and values of W
mydata = expand.grid(nsim=1:20,np=c(1000,5000,10000,20000),W=c(0.5,1,2),stringsAsFactors=FALSE)
require(doMC)
registerDoMC()
require(plyr)
rm.lmarglik = maply(.data = mydata, .fun = pf_lmarglik, .parallel = TRUE)
save(rm.lmarglik, file="../data/mc_1sim_test.rdata")