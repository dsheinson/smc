source("mc_functions.r")

# Simulate states and observations from dlm
F <- G <- V <- W <- sigma <- 1
nt = 100
sim = dlm.sim(nt, F, G, V, W, sigma)

# Calculate sufficient statistics of filtered distributions
post = dlm.post(sim$y, F, G, V, W, 1, 1, 0, 1)

# Plot true data, states, 95% CI of filtered states, and 95% one-step ahead PI
burn = 1
lk = qt(0.025,2*post$a)*sqrt(post$C[1,1,]*(post$b/post$a)) + post$m[,1]
uk = qt(0.975,2*post$a)*sqrt(post$C[1,1,]*(post$b/post$a)) + post$m[,1]
lm = qt(0.025,2*post$a[-nt])*sqrt(post$Q[1,1,]*(post$b[-nt]/post$a[-nt])) + post$f[,1]
um = qt(0.975,2*post$a[-nt])*sqrt(post$Q[1,1,]*(post$b[-nt]/post$a[-nt])) + post$f[,1]
gmin = min(sim$x, sim$y, lk[-(1:burn)], lm[-(1:burn)])
gmax = max(sim$x, sim$y, uk[-(1:burn)], um[-(1:burn)])

plot(0:nt,sim$x[,1],ylim=c(gmin,gmax),type="l",xlab=expression(t),ylab="Position")
points(1:nt,sim$y[,1])
lines(0:nt,lk,col=2)
lines(0:nt,uk,col=2)
lines(1:nt,lm,col=4)
lines(1:nt,um,col=4)
legend("bottomright",legend=c(expression(x,y),"95% CI","95% PI"),lty=c(1,NA,1,1),pch=c(NA,1,NA,NA),col=c(1,1,2,4))

# Calculate true log marginal likelihood of the data
dlm.lmarglik(sim$y[,1], post$f[,1], post$Q[1,1,], post$a, post$b)

# Run resample move particle filter and approximate log marginal likelihood of the data
source("rm_pf.r")
source("rm_test_functions.r")
np = 100
a0 = b0 = 1
rprior1 <- function(j) rprior(j,a0,b0)
rmove <- function(y, x, theta) rm_mcmc(y, x, theta, a0, b0, 1)
out = rm_pf(sim$y[,1], dllik, revo, rprior1, rmove, np, method="stratified", nonuniformity="ess", threshold=0.8, log=FALSE)
pf.lmarglik(out)
