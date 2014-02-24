source("dlm_sim.r")

# Set data path
dpath = "../data/"

# Simulate data
N = 20
nt = 100
F = G = V = W = C0 = 1
m0 = 0
mysims = list()
for(j in 1:N) mysims[[j]] = dlm.sim(nt, F, G, V, W, m0, C0)

# Save data
save(mysims, file=paste(dpath,"rw_sim.rdata",sep=""))
