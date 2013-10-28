source("mc_functions.r")

# Load simulated data
load("../data/dlm_sim.rdata")
nt = dim(sims$y)[2]

# Calculate true posterior distribution under correct model
post = dlm.post(sims$y[1,], sims[[3]]$F, sims[[3]]$G, sims[[3]]$V, sims[[3]]$W, 1, 1, 0, 1)

# Plot simulated data/states with 95\% CI and PI
burn = 1
lk = qt(0.025,2*post$a)*sqrt(post$C[1,1,]*(post$b/post$a)) + post$m[,1]
uk = qt(0.975,2*post$a)*sqrt(post$C[1,1,]*(post$b/post$a)) + post$m[,1]
lm = qt(0.025,2*post$a[-nt])*sqrt(post$Q[1,1,]*(post$b[-nt]/post$a[-nt])) + post$f[,1]
um = qt(0.975,2*post$a[-nt])*sqrt(post$Q[1,1,]*(post$b[-nt]/post$a[-nt])) + post$f[,1]
gmin = min(sims$x[1,], sims$y[1,], lk[-(1:burn)], lm[-(1:burn)])
gmax = max(sims$x[1,], sims$y[1,], uk[-(1:burn)], um[-(1:burn)])

pdf(file="../graphs/mc_1sim_test-states.pdf",width=10,height=5)
par(mfrow=c(1,2))
plot(0:nt,sims$x[1,],ylim=c(gmin,gmax),type="l",main="95% CI for Filtered States",xlab=expression(t),ylab="Position")
lines(0:nt,lk,col=2)
lines(0:nt,uk,col=2)
legend("bottomright",legend=c(expression(x),"95% CI"),lty=c(1,1),col=c(1,2))
plot(1:nt,sims$y[1,],ylim=c(gmin,gmax),type="l",main="95% One-Step Ahead PI",xlab=expression(t),ylab="")
lines(1:nt,lm,col=4)
lines(1:nt,um,col=4)
legend("bottomright",legend=c(expression(y),"95% PI"),lty=c(1,1),col=c(1,4))
dev.off()

# Plot 95% CI for filtered precision
lp = qgamma(0.025,post$a,post$b)
up = qgamma(0.975,post$a,post$b)
pdf(file="../graphs/mc_1sim_test-precision.pdf")
plot(0:nt,lp,type="l",col=2,ylim=c(min(lp,1/sims[[3]]$sigma^2),max(up,1/sims[[3]]$sigma^2)),main="95% CI for Filtered Precision",xlab=expression(t),ylab=expression(phi))
lines(0:nt,up,col=2)
abline(h=1/sims[[3]]$sigma^2)
legend("topright",legend=c("Truth","95% CI"),lty=c(1,1),col=c(1,2))
dev.off()

# Calculate true log marginal likelihoods under each model
true_lmarglik = function(W)
{
  post = dlm.post(sims$y[1,], sims[[3]]$F, sims[[3]]$G, sims[[3]]$V, W, 1, 1, 0, 1)
  return(dlm.lmarglik(sims$y[1,], post$f[,1], post$Q[1,1,], post$a, post$b))
}
require(plyr)
true.lmarglik = maply(data.frame(W=c(0.5,1,2)), true_lmarglik)

# Load approximate log marginal likelihoods for particle filter runs
load("../data/mc_1sim_test.rdata")
nsim = as.numeric(dimnames(rm.lmarglik)[[1]])
np = as.numeric(dimnames(rm.lmarglik)[[2]])
W = as.numeric(dimnames(rm.lmarglik)[[3]])

# Plot histograms of of log marginal likelihoods for pfs compared with truth
require(plyr)
max.densities = aaply(rm.lmarglik, 2:3, function(x) max(hist(x,plot=FALSE)$intensities,na.rm=TRUE))
min.densities = aaply(rm.lmarglik, 2:3, function(x) min(hist(x,plot=FALSE)$intensities,na.rm=TRUE))
max.breaks = aaply(rm.lmarglik, 2:3, function(x) max(hist(x,plot=FALSE)$breaks,na.rm=TRUE))
min.breaks = aaply(rm.lmarglik, 2:3, function(x) min(hist(x,plot=FALSE)$breaks,na.rm=TRUE))
ylabs = c("Density","","")
xlabs = c("Log Marginal Likelihood","","")
cols = gray.colors(length(np))
mains = expression(paste(tilde(W)," = 0.5",sep=""),paste(tilde(W)," = 1",sep=""),paste(tilde(W)," = 2",sep=""))
pdf(file="../graphs/mc_1sim_test-histograms.pdf",width=15,height=5)
par(mfrow=c(1,3),mar=c(5,6,4,2)+0.1)
for(i in 1:3)
{
  hist(rm.lmarglik[,1,i],breaks=4,freq=FALSE,col=cols[1],main=mains[i],xlab=xlabs[i],ylab=ylabs[i],xlim=c(min(min.breaks),max(max.breaks)),ylim=c(min(min.densities[,i]),max(max.densities[,i])),cex.axis=1.5,cex.lab=1.75,cex.main=2)
  hist(rm.lmarglik[,2,i],breaks=4,freq=FALSE,add=TRUE,col=cols[2])
  hist(rm.lmarglik[,3,i],breaks=4,freq=FALSE,add=TRUE,col=cols[3])
  hist(rm.lmarglik[,4,i],breaks=4,freq=FALSE,add=TRUE,col=cols[4])
  abline(v=true.lmarglik[i],lwd=2,col=2)
  mtext(round(true.lmarglik[i],2),side=1,at=true.lmarglik[i],cex=0.75)
  if(i == 1) legend("topleft",c(paste(np," particles",sep=""),"Truth"),fill=c(gray.colors(length(np)),NA),border=c(rep("black",length(np)),"white"),lty=c(rep(NA,length(np)),1),lwd=c(rep(NA,length(np)),2),col=c(rep(NA,length(np)),2))
}
dev.off()

# Compute true posterior model probabilities
true.postModProbs = postModProbs(true.lmarglik, rep(1/length(W), length(W)))

# Compute rm_pf approximate posterior model probabilities
rm.postModProbs = aaply(rm.lmarglik, 1:2, function(x) postModProbs(x, rep(1/length(W), length(W))))

# How often does the pf pick the "right" model?
rm.max = apply(rm.postModProbs, 1:2, function(x) which(x == max(x)))
prob.matchTruth = apply(rm.max, 2, function(x) sum(x == which(true.postModProbs == max(true.postModProbs)))/length(nsim))

# Ternary diagrams of posterior model probabilities
require(compositions)
pdf(file="../graphs/mc_1sim_test-ternary.pdf",width=10,height=10)
par(mfrow=c(2,2))
for(i in 1:length(np))
{
  plot(acomp(rm.postModProbs[,i,]),labels=c(expression(paste(tilde(W)," = 0.5",sep="")),expression(paste(tilde(W)," = 1",sep="")),expression(paste(tilde(W)," = 2",sep=""))),lwd=2)
  plot(acomp(true.postModProbs),add=TRUE,col=2,lwd=2)
  mtext(paste(np[i]," particles",sep=""),side=3,cex=2)
  if(i==1) legend("topleft",legend="True Posterior",pch=1,col=2,pt.lwd=2)
}
dev.off()