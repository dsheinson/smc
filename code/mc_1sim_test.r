source("mc_functions.r")
source("rm_pf.r")
source("rm_test_functions.r")

# Simulate univariate time series of states and observations
F <- G <- V <- W <- sigma <- 1
nt = 100
sim = dlm.sim(nt, F, G, V, W, sigma)

# Save simulated data
save.image("../data/mc_1sim_test-truth.rdata")

# Function to run resample move particle filter and approximate log marginal likelihood of the data
pf_lmarglik <- function(np, W, nrep, corr = TRUE)
{
  a0 = b0 = 1
  rprior1 <- function() rprior(a0,b0)
  revo1 <- function(x, theta) revo(x, theta, w = W)
  rmove <- function(y, x, theta) rm_mcmc(y, x, theta, a0, b0, 1)
  if(corr) set.seed(W)
  out = rm_pf(sim$y[,1], dllik, revo1, rprior1, rmove, np, progress = FALSE, method="stratified", nonuniformity="ess", threshold=0.8, log=FALSE)
  save(out, file=paste("../data/mc_1sim_test-",np,"-",nrep,"-",W,".rdata",sep=""))
  cat(np,nrep,W,'\n')
}

# Apply pf_lmarglik for many pf runs with different number of particles
mydata = expand.grid(np = c(1000, 5000, 10000, 20000), W = c(0.5, 1, 2), nrep = seq(1,20,1), stringsAsFactors=FALSE)
require(doMC)
registerDoMC()
require(plyr)
m_ply(.data = mydata, .fun = pf_lmarglik, .parallel = TRUE)