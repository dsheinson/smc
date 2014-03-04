source("dlm_sim.r")

# Set data path
dpath = "../data/"

# Simulate data
N = 20
nt = 200
F = G = V = W = C0 = 1
x0 = 0
mysims = list()
for(j in 1:N) mysims[[j]] = dlm.sim(nt, F, G, V, W, x0)

# Save data
save(mysims, file=paste(dpath,"rw_sim.rdata",sep=""))
