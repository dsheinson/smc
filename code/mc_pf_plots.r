source("mc_functions.r")
require(plyr)

# Load simulated data sets
load("../data/mc_pf_test-sims.rdata")

# load log marginal likelihoods
load("../data/mc_pf_lmarglik.rdata")
nsim = as.numeric(rownames(lmarglik.out$true.lmarglik))
W = as.numeric(colnames(lmarglik.out$true.lmarglik))
np = as.numeric(rownames(lmarglik.out$rm.lmarglik))
N = length(nsim)

# Compute true posterior model probabilities
true.postModProbs = apply(lmarglik.out$true.lmarglik, 1, function(x) postModProbs(x, rep(1/length(W), length(W))))

# Compute rm_pf approximate posterior model probabilities
rm.postModProbs = aaply(lmarglik.out$rm.lmarglik, 1:2, function(x) postModProbs(x, rep(1/length(W), length(W))))

# 3-way contingency tables
true.max = as.factor(apply(true.postModProbs, 2, function(x) which(x == max(x))))
levels(true.max) = seq(1,length(W),1)
rm.max = aaply(rm.postModProbs, 1:2, function(x) which(x == max(x)))
for(i in 1:dim(rm.max)[1])
{
  rm.max.f = as.factor(rm.max[i,])
  levels(rm.max.f) <- seq(1,length(W),1)
  tab = table(true.max, rm.max.f)
  dimnames(tab) = list(truth=c("M_1","M_2","M_3"),pf=c("M_1","M_2","M_3"))
  require(xtable)
  print(xtable(tab,caption=paste("Contingency table of models with highest posterior model probability among three models calculated analytically and approximated using the resample-move particle filter with ",np[i]," particles.",sep=""),label=paste("tab:cont-",np[i],sep="")),file=paste("../data/mc_pf_test-contTable-",np[i],".txt",sep=""),type="latex") 
}

# Plot ternary diagrams for posterior model probabilities
require(compositions)
grad = order(true.postModProbs[2,],decreasing=TRUE)
pdf(file="../graphs/mc_pf_test-ternary.pdf",width=10,height=10)
par(mfrow=c(2,2))
plot(acomp(t(true.postModProbs)[grad,]),labels=c(expression(paste(tilde(W)," = 0.5",sep="")),expression(paste(tilde(W)," = 1",sep="")),expression(paste(tilde(W)," = 2",sep=""))),col=gray.colors(N),lwd=2)
mtext("True Posterior",side=3,cex=2)
for(i in 1:length(np))
{
  plot(acomp(rm.postModProbs[i,grad,]),labels=c(expression(paste(tilde(W)," = 0.5",sep="")),expression(paste(tilde(W)," = 1",sep="")),expression(paste(tilde(W)," = 2",sep=""))),col=gray.colors(N),lwd=2)
  mtext(paste(np[i]," particles",sep=""),side=3,cex=2)
}
dev.off()

# Plot binary diagrams for posterior model probabilities
# Comparing W = 0.5 versus W = 1.0
pdf(file="../graphs/mc_pf_test-binary.pdf",width=12,height=8)
par(mfrow=c(2,length(np)),mar=c(5,6,4,1)+.1)
ylabs=c(expression(paste("P(",tilde(W)," = 1) vs P(",tilde(W)," = 0.5)",sep="")),rep("",length(np)-1))
true.postModProbs = apply(lmarglik.out$true.lmarglik[,1:2], 1, function(x) postModProbs(x, rep(1/length(W[1:2]), length(W[1:2]))))
rm.postModProbs = aaply(lmarglik.out$rm.lmarglik[,,1:2], 1:2, function(x) postModProbs(x, rep(1/length(W[1:2]), length(W[1:2]))))
grad = order(true.postModProbs[2,],decreasing=TRUE)
for(i in 1:length(np))
{
  plot(rep(1,N),true.postModProbs[2,grad],main=paste(np[i]," particles",sep=""),col=gray.colors(N),xaxt = 'n',yaxt='n',lwd=2,xlim=c(0.9,2.1),ylim=c(0,1),xlab="",ylab=ylabs[i],cex.lab=1.75,cex.main=1.75)
  axis(side=1, at=c(1,2), labels=c("Direct Calculation","PF approx."),cex.axis=1.65)
  axis(side=2, at=c(0,0.5,1),cex.axis=1.65)
  points(rep(2,N),rm.postModProbs[i,grad,2],col=gray.colors(N),lwd=2)
  segments(rep(1,N), true.postModProbs[2,grad], rep(2,N), rm.postModProbs[i,grad,2],col=gray.colors(N))
}
ylabs=c(expression(paste("P(",tilde(W)," = 1) vs P(",tilde(W)," = 2)",sep="")),rep("",length(np)-1))
true.postModProbs = apply(lmarglik.out$true.lmarglik[,2:3], 1, function(x) postModProbs(x, rep(1/length(W[2:3]), length(W[2:3]))))
rm.postModProbs = aaply(lmarglik.out$rm.lmarglik[,,2:3], 1:2, function(x) postModProbs(x, rep(1/length(W[2:3]), length(W[2:3]))))
grad = order(true.postModProbs[1,],decreasing=TRUE)
for(i in 1:length(np))
{
  plot(rep(1,N),true.postModProbs[1,grad],col=gray.colors(N),xaxt = 'n',yaxt='n',lwd=2,xlim=c(0.9,2.1),ylim=c(0,1),xlab="",ylab=ylabs[i],cex.lab=1.75,cex.main=1.75)
  axis(side=1, at=c(1,2), labels=c("Direct Calculation","PF approx."),cex.axis=1.65)
  axis(side=2, at=c(0,0.5,1),cex.axis=1.65)
  points(rep(2,N),rm.postModProbs[i,grad,1],col=gray.colors(N),lwd=2)
  segments(rep(1,N), true.postModProbs[1,grad], rep(2,N), rm.postModProbs[i,grad,1],col=gray.colors(N))
}
dev.off()