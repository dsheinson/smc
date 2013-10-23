source("rm_pf.r")
source("rm_test_functions.r")

load("../data/mc_pf_test-sims.rdata")

mc_pf_test <- function(np, nsim, W)
{
  revo1 <- function(x, theta) rnorm(1,x,sqrt(theta*W))
  a0 = b0 = 1
  rprior1 <- function(j) rprior(j,a0,b0)
  rmove <- function(y, x, theta) rm_mcmc(y, x, theta, a0, b0, 1)
  out = rm_pf(sims$y[nsim,], dllik, revo1, rprior1, rmove, np, progress = FALSE, method="stratified", nonuniformity="ess", threshold=0.8, log=FALSE)
  save(out, file=paste("../data/mc_pf_test-",np,"-",nsim,"-",W,".rdata",sep=""))
  cat(np,nsim,W,'\n')
}

mydata = expand.grid(np=c(1000, 5000,10000),nsim=seq(1,20,1),W=c(.5,1,2),stringsAsFactors=FALSE)
require(doMC)
registerDoMC()
require(plyr)
m_ply(.data = mydata, .fun = mc_pf_test, .parallel = TRUE)