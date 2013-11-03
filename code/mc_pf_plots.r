source("mc_functions.r")
require(plyr)

# Load simulated data and approximate log marginal likelihoods
load("../data/dlm_sim.rdata")
nt = dim(sims$y)[2]
load("../data/mc_pf_test.rdata")
nsim = as.numeric(dimnames(rm.lmarglik)[[1]])
np = as.numeric(dimnames(rm.lmarglik)[[2]])
W = as.numeric(dimnames(rm.lmarglik)[[3]])
N = length(nsim)

# Calculate true log marginal likelihoods under each model
true_lmarglik = function(nsim, W)
{
  post = dlm.post(sims$y[nsim,], sims[[3]]$F, sims[[3]]$G, sims[[3]]$V, W, 1, 1, 0, 1)
  return(dlm.lmarglik(sims$y[nsim,], post$f[,1], post$Q[1,1,], post$a, post$b))
}
true.lmarglik = maply(expand.grid(nsim=nsim, W=W, stringsAsFactors=TRUE), true_lmarglik)

# Compute true posterior model probabilities
true.postModProbs = aaply(true.lmarglik, 1, function(x) postModProbs(x, rep(1/length(W), length(W))))

# Compute rm_pf approximate posterior model probabilities
rm.postModProbs = aaply(rm.lmarglik, 1:2, function(x) postModProbs(x, rep(1/length(W), length(W))))

# 3-way contingency tables
true.max = factor(apply(true.postModProbs, 1, function(x) which(x == max(x))),levels=seq(1,length(W),1))
rm.max = aaply(rm.postModProbs, 1:2, function(x) which(x == max(x)))
for(i in 1:dim(rm.max)[2])
{
  rm.max.f = factor(rm.max[,i],levels=seq(1,length(W),1))
  tab = table(true.max, rm.max.f)
  dimnames(tab) = list(truth=c("M_1","M_2","M_3"),pf=c("M_1","M_2","M_3"))
  require(xtable)
  print(xtable(tab,caption=paste("Contingency table of models with highest posterior model probability among three models calculated analytically and approximated using the resample-move particle filter with ",np[i]," particles.",sep=""),label=paste("tab:cont-",np[i],sep="")),file=paste("../data/mc_pf_test-contTable-",np[i],".txt",sep=""),type="latex") 
}

# Plot ternary diagrams for posterior model probabilities
require(compositions)
grad = order(true.postModProbs[,2],decreasing=TRUE)
labels = rep(NA,length(W))
for(i in 1:length(W)) labels[i] = as.expression(bquote(paste(tilde(W)," = ",sep="") ~ .(W[i])))
pdf(file="../graphs/mc_pf_test-ternary.pdf",width=10,height=10)
par(mfrow=c(2,2))
plot(acomp(true.postModProbs[grad,]),labels=labels,col=gray.colors(N),lwd=2)
mtext("True Posterior",side=3,cex=2)
for(i in 1:length(np))
{
  plot(acomp(rm.postModProbs[grad,i,]),labels=labels,col=gray.colors(N),lwd=2)
  mtext(paste(np[i]," particles",sep=""),side=3,cex=2)
}
dev.off()

# Plot binary diagrams for posterior model probabilities
pdf(file="../graphs/mc_pf_test-binary.pdf",width=12,height=8)
par(mfrow=c(2,length(np)),mar=c(5,6,4,1)+.1)
ylabs=c(substitute(paste("P(",tilde(W)," = ",aa,") vs P(",tilde(W)," = ",ab,")",sep=""),list(aa=W[2],ab=W[1])),rep("",length(np)-1))
true.postModProbs = aaply(true.lmarglik[,1:2], 1, function(x) postModProbs(x, rep(1/length(W[1:2]), length(W[1:2]))))
rm.postModProbs = aaply(rm.lmarglik[,,1:2], 1:2, function(x) postModProbs(x, rep(1/length(W[1:2]), length(W[1:2]))))
grad = order(true.postModProbs[,2],decreasing=TRUE)
for(i in 1:length(np))
{
  plot(rep(1,N),true.postModProbs[grad,2],main=paste(np[i]," particles",sep=""),col=gray.colors(N),xaxt = 'n',yaxt='n',lwd=2,xlim=c(0.9,2.1),ylim=c(0,1),xlab="",ylab=ylabs[i],cex.lab=1.75,cex.main=1.75)
  axis(side=1, at=c(1,2), labels=c("Direct Calculation","PF approx."),cex.axis=1.65)
  axis(side=2, at=c(0,0.5,1),cex.axis=1.65)
  points(rep(2,N),rm.postModProbs[grad,i,2],col=gray.colors(N),lwd=2)
  segments(rep(1,N), true.postModProbs[grad,2], rep(2,N), rm.postModProbs[grad,i,2],col=gray.colors(N))
}
ylabs=c(substitute(paste("P(",tilde(W)," = ",aa,") vs P(",tilde(W)," = ",ab,")",sep=""),list(aa=W[2],ab=W[3])),rep("",length(np)-1))
true.postModProbs = aaply(true.lmarglik[,2:3], 1, function(x) postModProbs(x, rep(1/length(W[2:3]), length(W[2:3]))))
rm.postModProbs = aaply(rm.lmarglik[,,2:3], 1:2, function(x) postModProbs(x, rep(1/length(W[2:3]), length(W[2:3]))))
grad = order(true.postModProbs[,1],decreasing=TRUE)
for(i in 1:length(np))
{
  plot(rep(1,N),true.postModProbs[grad,1],col=gray.colors(N),xaxt = 'n',yaxt='n',lwd=2,xlim=c(0.9,2.1),ylim=c(0,1),xlab="",ylab=ylabs[i],cex.lab=1.75,cex.main=1.75)
  axis(side=1, at=c(1,2), labels=c("Direct Calculation","PF approx."),cex.axis=1.65)
  axis(side=2, at=c(0,0.5,1),cex.axis=1.65)
  points(rep(2,N),rm.postModProbs[grad,i,1],col=gray.colors(N),lwd=2)
  segments(rep(1,N), true.postModProbs[grad,1], rep(2,N), rm.postModProbs[grad,i,1],col=gray.colors(N))
}
dev.off()