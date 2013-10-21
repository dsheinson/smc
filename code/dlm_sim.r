source("mc_functions.r")

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