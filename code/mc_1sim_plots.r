# Load simulated data
load("../data/mc_1sim_test-truth.rdata")

# Calculate true log marginal likelihood of the simulated data
source("mc_functions.r")
post = dlm.post(sim$y, F, G, V, W, 1, 1, 0, 1)
true.marg = dlm.lmarglik(sim$y[,1], post$f[,1], post$Q[1,1,], post$a, post$b)

# Plot simulated data/states with 95\% CI and PI
burn = 1
lk = qt(0.025,2*post$a)*sqrt(post$C[1,1,]*(post$b/post$a)) + post$m[,1]
uk = qt(0.975,2*post$a)*sqrt(post$C[1,1,]*(post$b/post$a)) + post$m[,1]
lm = qt(0.025,2*post$a[-nt])*sqrt(post$Q[1,1,]*(post$b[-nt]/post$a[-nt])) + post$f[,1]
um = qt(0.975,2*post$a[-nt])*sqrt(post$Q[1,1,]*(post$b[-nt]/post$a[-nt])) + post$f[,1]
gmin = min(sim$x, sim$y, lk[-(1:burn)], lm[-(1:burn)])
gmax = max(sim$x, sim$y, uk[-(1:burn)], um[-(1:burn)])

pdf(file="../graphs/mc_1sim_test-states.pdf",width=10,height=5)
par(mfrow=c(1,2))
plot(0:nt,sim$x[,1],ylim=c(gmin,gmax),type="l",main="95% CI for Filtered States",xlab=expression(t),ylab="Position")
lines(0:nt,lk,col=2)
lines(0:nt,uk,col=2)
legend("bottomright",legend=c(expression(x),"95% CI"),lty=c(1,1),col=c(1,2))
plot(1:nt,sim$y[,1],ylim=c(gmin,gmax),type="l",main="95% One-Step Ahead PI",xlab=expression(t),ylab="")
lines(1:nt,lm,col=4)
lines(1:nt,um,col=4)
legend("bottomright",legend=c(expression(y),"95% PI"),lty=c(1,1),col=c(1,4))
dev.off()

# Plot 95% CI for filtered precision
lp = qgamma(0.025,post$a,post$b)
up = qgamma(0.975,post$a,post$b)
pdf(file="../graphs/mc_1sim_test-precision.pdf")
plot(0:nt,lp,type="l",col=2,ylim=c(min(lp,1/sigma^2),max(up,1/sigma^2)),main="95% CI for Filtered Precision",xlab=expression(t),ylab=expression(phi))
lines(0:nt,up,col=2)
abline(h=1/sigma^2)
legend("topright",legend=c("Truth","95% CI"),lty=c(1,1),col=c(1,2))
dev.off()

# Load approximate log marginal likelihoods for particle filter runs
np = c(100, 500, 1000)
N = 20
pf.margs = matrix(NA, nr=N, nc=length(np))
for(i in 1:N)
{
  for(j in 1:length(np))
  {
    load(paste("../data/mc_1sim_test-",np[j],"-",20*(j-1)+i,".rdata",sep=""))
    pf.marg = pf.lmarglik(out)
    pf.margs[i,j] = pf.marg
  }
}

# Plot histograms of of log marginal likelihoods for pfs compared with truth
pdf(file="../graphs/mc_1sim_hist.pdf")
hist(pf.margs[,1],breaks=5,main="Histogram of Log Marginal Likelihoods",xlab="Log Marginal Likelihood",xlim=c(min(pf.margs,true.marg),max(pf.margs,true.marg)))
hist(pf.margs[,2],breaks=5,add=TRUE,col="gray")
hist(pf.margs[,3],breaks=5,add=TRUE,col="black")
abline(v=true.marg,lwd=2,col=2)
legend(-370,10,c("100 particles","500 particles","1000 particles","Truth"),fill=c("white","gray","black",NA),border=c("black","black","black","white"),lty=c(NA,NA,NA,1),lwd=c(NA,NA,NA,2),col=c(NA,NA,NA,2),cex=0.85)
dev.off()