source("mc_functions.r")
source("rm_pf.r")
source("rm_test_functions.r")

# Simulate data
N = 20
nt = 100
F = G = V = W = sigma = 1
x = matrix(NA, nr = N, nc = nt + 1)
y = matrix(NA, nr = N, nc = nt)
for(j in 1:N)
{
  sim = dlm.sim(nt, F, G, V, W, sigma)
  x[j,] = sim$x[,1]
  y[j,] = sim$y[,1]
}
sims = list(x=x,y=y)
save(sims, file="../data/mc_pf_test-sims.rdata")

mc_pf_test <- function(np, nsim, W)
{
  revo1 <- function(x, theta) rnorm(1,x,sqrt(theta*W))
  a0 = b0 = 1
  rprior1 <- function(j) rprior(j,a0,b0)
  rmove <- function(y, x, theta) rm_mcmc(y, x, theta, a0, b0, 1)
  out = rm_pf(y[nsim,], dllik, revo1, rprior1, rmove, np, progress = FALSE, method="stratified", nonuniformity="ess", threshold=0.8, log=FALSE)
  save(out, file=paste("../data/mc_pf_test-",np,"-",nsim,"-",W,".rdata",sep=""))
  print(paste(np,nsim,W,sep="-"))
}

require(plyr)
mydata = expand.grid(np=c(5,10,15),nsim=seq(1,5,1),W=c(.5,1,2),stringsAsFactors=TRUE)
m_ply(mydata, mc_pf_test)
#mydata = expand.grid(np=c(100,500,1000),nsim=seq(1,20,1),W=c(.5,1,2),stringsAsFactors=FALSE)
#require(doMC)
#registerDoMC()
#m_ply(mydata, mc_pf_test, .parallel=TRUE)