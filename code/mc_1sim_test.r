source("mc_functions.r")
source("rm_pf.r")
source("rm_test_functions.r")

# Simulate univariate time series of states and observations
F <- G <- V <- W <- sigma <- 1
nt = 100
sim = dr.sim(nt, F, G, V, W, sigma)

# Save simulated data
save.image("../data/mc_1sim_test-truth.rdata")

# Function to run resample move particle filter and approximate log marginal likelihood of the data
pf_lmarglik <- function(np, label)
{
  a0 = b0 = 1
  rprior1 <- function(j) rprior(j,a0,b0)
  rmove <- function(y, x, theta) rm_mcmc(y, x, theta, a0, b0, 1)
  out = rm_pf(sim$y[,1], dllik, revo, rprior1, rmove, np, progress = FALSE, method="stratified", nonuniformity="ess", threshold=0.8, log=FALSE)
  print(np)
  pf.marg = pf.lmarglik(out)
  save(pf.marg, file=paste("../data/mc_1sim_test-",np,"-",label,".rdata",sep=""))
}

# Apply pf_lmarglik for many pf runs with different number of particles
mydata = data.frame(np = rep(c(100, 500, 1000), rep(20, 3)), label=seq(1,60,1))
require(plyr)
m_ply(.data = mydata, .fun = pf_lmarglik)