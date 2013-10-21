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
np = c(5,10,15)
rm_lmarglik = function(np, nsim, W)
{
  load(paste("../data/mc_pf_test-",np,"-",nsim,"-",W,".rdata",sep=""))
  return(pf.lmarglik(out))
}
rm.lmarglik = maply(expand.grid(np=np,nsim=seq(1,5,1),W=W.test), rm_lmarglik)

# Compute true posterior model probabilities
true.postModProbs = apply(true.lmarglik, 1, function(x) postModProbs(x, rep(1/length(W.test), length(W.test))))

# Compute rm_pf approximate posterior model probabilities
rm.postModProbs = aaply(rm.lmarglik, 1:2, function(x) postModProbs(x, rep(1/length(W.test), length(W.test))))

# Contingency tables
r = 1:5
true.max = apply(true.postModProbs[,r], 2, function(x) which(x == max(x)))
rm.max = aaply(rm.postModProbs[,r,], 1:2, function(x) which(x == max(x)))
for(i in 1:dim(rm.max)[1])
{
  tab = table(true.max, rm.max[i,])
#  dimnames(tab) = list(Truth=c("M_1","M_2","M_3"),Approx=c("M_1","M_2","M_3"))
  require(xtable)
  print(xtable(tab),file=paste("../data/mc_pf_test-contTable-",np[i],".txt",sep=""),type="latex") 
}

# Plot ternary diagrams for posterior model probabilities
require(compositions)
grad = order(true.postModProbs[2,],decreasing=TRUE)
pdf(file="../graphs/10-23-13/mc_pf_test-ternary.pdf",width=10,height=10)
par(mfrow=c(2,2))
plot(acomp(t(true.postModProbs)[grad,]),labels=c("W = 0.5","W = 1","W = 2"),col=gray.colors(N),lwd=2)
mtext("True Posterior",side=3,cex=2)
for(i in 1:length(np))
{
  plot(acomp(rm.postModProbs[i,r,]),labels=c("W = 0.5","W = 1","W = 2"),col=gray.colors(N)[grad[r]],lwd=2)
  mtext(paste(np[i]," particles",sep=""),side=3,cex=2)
}
dev.off()

# Plot binary diagrams for posterior model probabilities
# Comparing W = 0.5 versus W = 1.0
pdf(file="../graphs/10-23-13/mc_pf_test-binary.pdf",width=12,height=8)
par(mfrow=c(2,length(np)),mar=c(5,6,4,1)+.1)
ylabs=c(expression(paste("P(",tilde(W)," = 1) vs P(",tilde(W)," = 0.5)",sep="")),rep("",length(np)-1))
true.postModProbs = apply(true.lmarglik[,1:2], 1, function(x) postModProbs(x, rep(1/length(W.test[1:2]), length(W.test[1:2]))))
rm.postModProbs = aaply(rm.lmarglik[,,1:2], 1:2, function(x) postModProbs(x, rep(1/length(W.test[1:2]), length(W.test[1:2]))))
grad = order(true.postModProbs[2,],decreasing=TRUE)
for(i in 1:length(np))
{
  plot(rep(1,N),true.postModProbs[2,grad],main=paste(np[i]," particles",sep=""),col=gray.colors(N),xaxt = 'n',yaxt='n',lwd=2,xlim=c(0.9,2.1),ylim=c(0,1),xlab="",ylab=ylabs[i],cex.lab=1.75,cex.main=1.75)
  axis(side=1, at=c(1,2), labels=c("Direct Calculation","PF approx."),cex.axis=1.65)
  axis(side=2, at=c(0,0.5,1),cex.axis=1.65)
  points(rep(2,5),rm.postModProbs[i,r,2],col=gray.colors(N)[grad[r]],lwd=2)
  segments(rep(1,5), true.postModProbs[2,r], rep(2,5), rm.postModProbs[i,r,2])
}
ylabs=c(expression(paste("P(",tilde(W)," = 1) vs P(",tilde(W)," = 2)",sep="")),rep("",length(np)-1))
true.postModProbs = apply(true.lmarglik[,2:3], 1, function(x) postModProbs(x, rep(1/length(W.test[2:3]), length(W.test[2:3]))))
rm.postModProbs = aaply(rm.lmarglik[,,2:3], 1:2, function(x) postModProbs(x, rep(1/length(W.test[2:3]), length(W.test[2:3]))))
grad = order(true.postModProbs[1,],decreasing=TRUE)
for(i in 1:length(np))
{
  plot(rep(1,N),true.postModProbs[1,grad],col=gray.colors(N),xaxt = 'n',yaxt='n',lwd=2,xlim=c(0.9,2.1),ylim=c(0,1),xlab="",ylab=ylabs[i],cex.lab=1.75,cex.main=1.75)
  axis(side=1, at=c(1,2), labels=c("Direct Calculation","PF approx."),cex.axis=1.65)
  axis(side=2, at=c(0,0.5,1),cex.axis=1.65)
  points(rep(2,5),rm.postModProbs[i,r,1],col=gray.colors(N)[grad[r]],lwd=2)
  segments(rep(1,5), true.postModProbs[1,r], rep(2,5), rm.postModProbs[i,r,1])
}
dev.off()