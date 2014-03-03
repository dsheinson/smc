source("fmri_mcmc_functions.r")

# Set data and graphics path
dpath = "../data/"
gpath = "../graphs/"

# Load data
load(paste(dpath,"dlm_ar_sim-20-M101.rdata",sep=""))
y = mysims[[1]]$y

# Set known values
V = 1
U = mysims[[1]]$true.params$U
F = mysims[[1]]$true.params$F
psi = list(V=V,U=U,F=F)

# Set initial values to truth and set prior on phi
initial = list(x=mysims[[1]]$x, theta = list(beta = mysims[[1]]$true.params$beta, sigma2m = mysims[[1]]$true.params$V[1,1], phi = mysims[[1]]$true.params$G[1,1], sigma2s = mysims[[1]]$true.params$W[1,1]))
prior = list(m0 = 0, b0 = c(0,0), B0 = 1e6*diag(2), am0 = 1e-6, bm0 = 1e-6, phi0 = 0, Phi0 = 1e6, am0 = 1e-6, bm0 = 1e-6)

mcmc.details = list(n.sims = 1100, n.thin = 1, n.burn = 100)
out = fmri_mcmc(y, psi, prior, initial, mcmc.details)