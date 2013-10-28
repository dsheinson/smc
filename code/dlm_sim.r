source("mc_functions.r")

# Simulate data
N = 20
nt = 100
F = G = V = W = sigma = 1
true.params = list(F=F,G=G,V=V,W=W,sigma=sigma)
x = matrix(NA, nr = N, nc = nt + 1)
y = matrix(NA, nr = N, nc = nt)
for(j in 1:N)
{
  sim = dlm.sim(nt, F, G, V, W, sigma)
  x[j,] = sim$x[,1]
  y[j,] = sim$y[,1]
}

# Save data
sims = list(x=x,y=y,true.params=true.params)
save(sims, file=paste("../data/dlm_sim.rdata",sep=""))