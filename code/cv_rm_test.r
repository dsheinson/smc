source("dlm_cv_functions.r")
source("pf_mc_functions.r")

# Set data and graphics paths
dpath = "../data/"
gpath = "../graphs/"

# Load data
load(paste(dpath,"rw_sim.rdata",sep=""))

# Function to test cv filter and log marginal likelihood calculation
cv_test <- function(n.sim, alpha = 0.05, burn = 1)
{
  F = mysims[[n.sim]]$true.params$F
  G = mysims[[n.sim]]$true.params$G
  V = mysims[[n.sim]]$true.params$V
  W = mysims[[n.sim]]$true.params$W
  a0 = 1
  b0 = 1
  m0 = 0
  C0 = 1
  nt = dim(mysims[[n.sim]]$y)[2]

  # Calculate sufficient statistics of filtered distributions
  post = cv.post(mysims[[n.sim]]$y, F, G, V, W, a0, b0, m0, C0)

  # Calculatue 95% credible intervals for filtered states, precision, and one-step ahead predictions
  lk = qt(alpha/2,2*post$a)*sqrt(post$C[1,1,]*(post$b/post$a)) + post$m[1,]
  uk = qt(1 - alpha/2,2*post$a)*sqrt(post$C[1,1,]*(post$b/post$a)) + post$m[1,]
  lm = qt(alpha/2,2*post$a[-nt])*sqrt(post$Q[1,1,]*(post$b[-nt]/post$a[-nt])) + post$f[1,]
  um = qt(1-alpha/2,2*post$a[-nt])*sqrt(post$Q[1,1,]*(post$b[-nt]/post$a[-nt])) + post$f[1,]
  lp = qgamma(alpha/2,post$a,post$b)
  up = qgamma(1-alpha/2,post$a,post$b)
  
  # Plot true data, states, 95% CI of filtered states, and 95% one-step ahead PI
  gmin = min(mysims[[n.sim]]$x, mysims[[n.sim]]$y, lk[-(1:burn)], lm[-(1:burn)])
  gmax = max(mysims[[n.sim]]$x, mysims[[n.sim]]$y, uk[-(1:burn)], um[-(1:burn)])
  pdf(paste(gpath,"cv-test-states-",n.sim,".pdf",sep=""))
  plot(0:nt,mysims[[n.sim]]$x[1,],ylim=c(gmin,gmax),type="l",xlab=expression(t),ylab="Position")
  points(1:nt,mysims[[n.sim]]$y[1,])
  lines(0:nt,lk,col=2)
  lines(0:nt,uk,col=2)
  lines(1:nt,lm,col=4)
  lines(1:nt,um,col=4)
  legend("bottomright",legend=c(expression(x,y),"95% CI","95% PI"),lty=c(1,NA,1,1),pch=c(NA,1,NA,NA),col=c(1,1,2,4))
  dev.off()
  
  # Plot 95% CI for filtered precision
  pdf(paste(gpath,"cv-test-precision-",n.sim,".pdf",sep=""))
  plot(0:nt,lp,type="l",col=2,ylim=c(min(lp,1),max(up,1)),main="95% CI for Filtered Precision",xlab=expression(t),ylab=expression(phi))
  lines(0:nt,up,col=2)
  abline(h=1)
  legend("topright",legend=c("Truth","95% CI"),lty=c(1,1),col=c(1,2))
  dev.off()
  
  # Check proportion of states outside 95% credible intervals
  pci = sum(uk < mysims[[n.sim]]$x[1,] | lk > mysims[[n.sim]]$x[1,]) / (nt+1)

  # Check proportion of data points outside 95% prediction intervals
  ppi = sum(um < mysims[[n.sim]]$y[1,] | lm > mysims[[n.sim]]$y[1,]) / nt

  # Calculate true log marginal likelihood of the data
  lmlik = cv.lmarglik(mysims[[n.sim]]$y, post$f, post$Q, post$a, post$b)
 
  return(list(pci=pci,ppi=ppi,lmlik=lmlik))
}

require(plyr)
cv.test = mlply(data.frame(n.sim=1:length(mysims)), cv_test)

# Graph percentage of points/states outside 95% C/P intervals
pc <- sapply(cv.test, function(x) x$pci)
pp <- sapply(cv.test, function(x) x$ppi)
N = length(mysims)
alpha = 0.05
se = sqrt(alpha*(1-alpha) / N)
gmax = max(alpha, pc+1.96*se, pp+1.96*se)+0.1
gmin = min(alpha, pc-1.96*se, pp-1.96*se)
pdf(paste(gpath,"cv-test-intervals.pdf",sep=""))
plot(1:N, pc, xlab = "Simulation", xaxt="n", ylab = "% outside interval", ylim = c(gmin,gmax), col=2)
axis(1, 1:N, 1:N)
points(1:N+.1, pp, col = 4)
segments(1:N, pc-1.96*se, 1:N, pc+1.96*se, col = 2)
segments(1:N+.1, pp-1.96*se, 1:N+.1, pp+1.96*se, col = 4)
abline(h=c(0,alpha), lty = c(1,2))
legend("topleft",c("CI","PI","nominal level"), lty = c(1,1,2), col=c(2,4,1))
title(paste(100*(1-alpha),"% CI for rejection rate"))
dev.off()

# Run resample move particle filter and approximate log marginal likelihood of the data
source("rm_pf.r")
source("rm_cv_functions.r")
source("pf_mc_functions.r")
n.sim = 1
F = mysims[[n.sim]]$true.params$F
G = mysims[[n.sim]]$true.params$G
V = mysims[[n.sim]]$true.params$V
W = mysims[[n.sim]]$true.params$W
m0 = 0
C0 = 1
a0 = b0 = 1
mydlm = list(F=F[1,1,1],G=G[1,1],V=V[1,1],W=W[1,1],m0=m0,C0=C0)
rmove <- function(y, x, theta) rm_mcmc(y, x, theta, a0, b0, mydlm, 1)
np = 100
out = rm_pf(mysims[[n.sim]]$y, dllik, revo, rprior, rmove, np, method="stratified", nonuniformity="ess", threshold=0.8, log=FALSE)
pf.lmarglik(out)

# Calculate 95% CI of states and precision
source("pf_functions.r")
state.quant = pf.quantile(out$state, out$weight, function(x, param=1) x, c(alpha/2,1-alpha/2))
theta.quant = pf.quantile(1/out$theta, out$weight, function(x, param=1) x, c(alpha/2,1-alpha/2))
# Compare with true posterior
post = cv.post(mysims[[n.sim]]$y, F, G, V, W, 1, 1, m0, C0)
lk = qt(alpha/2,2*post$a)*sqrt(post$C[1,1,]*(post$b/post$a)) + post$m[1,]
uk = qt(1 - alpha/2,2*post$a)*sqrt(post$C[1,1,]*(post$b/post$a)) + post$m[1,]
lp = qgamma(alpha/2,post$a,post$b)
up = qgamma(1-alpha/2,post$a,post$b)

# Plot 95% CI of filtered states
burn=1
nt = dim(out$state)[3] - 1
gmin = min(mysims[[n.sim]]$x, mysims[[n.sim]]$y, lk[-(1:burn)], state.quant[-(1:burn),1,1])
gmax = max(mysims[[n.sim]]$x, mysims[[n.sim]]$y, uk[-(1:burn)], state.quant[-(1:burn),1,2])
pdf(paste(gpath,"rm-test-states-",n.sim,".pdf",sep=""))
plot(0:nt,mysims[[n.sim]]$x[1,],ylim=c(gmin,gmax),type="l",xlab=expression(t),ylab="Position")
lines(0:nt,lk,col=2)
lines(0:nt,uk,col=2)
lines(0:nt,state.quant[,1,1],col=4)
lines(0:nt,state.quant[,1,2],col=4)
legend("bottomright",legend=c("PF Approx", "True Post", "True Sim"),lty=c(1,1,1),col=c(4,2,1))
title("95% CI for filtered states")
dev.off()

# Plot 95% CI for filtered precision
pdf(paste(gpath,"rm-test-precision-",n.sim,".pdf",sep=""))
plot(0:nt,lp,type="l",col=2,ylim=c(min(lp,1),max(up,1)),main="95% CI for Filtered Precision",xlab=expression(t),ylab=expression(phi))
lines(0:nt,up,col=2)
lines(0:nt,theta.quant[,1,1],col=4)
lines(0:nt,theta.quant[,1,2],col=4)
abline(h=1)
legend("topright",legend=c("PF Approx", "True Post", "True Sim"),lty=c(1,1,1),col=c(4,2,1))
dev.off()