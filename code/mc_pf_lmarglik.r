source("mc_functions.r")

# Load simulated data sets
load("../data/mc_pf_test-sims.rdata")

# Calculate true log marginal likelihoods for simulated data under different models
N = dim(sims$y)[1]
F = G = V = 1
W.test = c(0.5, 1, 2)
true_lmarglik = function(nsim, W)
{
  post = dlm.post(sims$y[nsim,], F, G, V, W, 1, 1, 0, 1)
  return(dlm.lmarglik(sims$y[nsim,], post$f[,1], post$Q[1,1,], post$a, post$b))
}
require(plyr)
true.lmarglik = maply(expand.grid(nsim=seq(1,N,1),W=W.test), true_lmarglik)

# Approximate (using resample-move particle filter) log marginal likelihoods for simulated data under different models
np = c(5000,10000,20000)
nsim = seq(1,20,1)
rm_lmarglik = function(np, nsim, W)
{
  file = paste("../data/mc_pf_test-",np,"-",nsim,"-",W,".rdata",sep="")
  load(file)
  print(file)
  return(pf.lmarglik(out))
}
rm.lmarglik = maply(expand.grid(np=np,nsim=nsim,W=W.test), rm_lmarglik)

# Save log marginal likelihoods
lmarglik.out = list(true.lmarglik=true.lmarglik,rm.lmarglik=rm.lmarglik)
save(lmarglik.out, file="../data/mc_pf_lmarglik.rdata")