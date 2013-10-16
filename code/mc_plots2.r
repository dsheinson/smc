# Plot histograms of of log marginal likelihoods for pfs
load("../data/mc_1pf_test-truth.rdata")
np = c(100, 500)
N = 20
pf.margs = matrix(NA, nr=N, nc=length(np))
for(i in 1:N)
{
  for(j in 1:length(np))
  {
    load(paste("../data/10-16-13/mc_1pf_test-",np[j],"-",20*(j-1)+i,".rdata",sep=""))
    pf.margs[i,j] = pf.marg
  }
}
hist(pf.margs[,1],breaks=5,main="Histogram of Log Marginal Likelihoods",xlab="Log Marginal Likelihood",xlim=c(min(pf.margs,true.marg),max(pf.margs,true.marg)))
abline(v=true.marg,lwd=2,col=2)
hist(pf.margs[,2],breaks=5,add=TRUE,col=8)
legend(-370,10,c("100 particles","500 particles","Truth"),fill=c("white","gray",NA),border=c("black","black","white"),lty=c(NA,NA,1),lwd=c(NA,NA,2),col=c(NA,NA,2),cex=0.85)

# Plot compositional data for posterior model probabilities
nt = 150
np = 100
N = 100
load(paste("../data/10-9-13/mc_pf_test-",nt,"-",np,"-",N,".rdata",sep=""))

# Comparing all three models
require(compositions)
windows(width=10,height=10)
par(mfcol=c(2,2))
plot(acomp(out$mod.post[,,1]),labels=c("W = 0.5","W = 1","W = 1.5"),col=rainbow(N),lwd=2)
mtext("True Posterior",side=2,line=4)
mtext(paste(np," particles",sep=""),side=3)
plot(acomp(out$mod.post[,,2]),labels=c("W = 0.5","W = 1","W = 1.5"),col=rainbow(N),lwd=2)
mtext("PF Approximation",side=2,line=4)
np = 200
load(paste("../data/10-9-13/mc_pf_test-",nt,"-",np,"-",N,".rdata",sep=""))
plot(acomp(out$mod.post[,,1]),labels=c("W = 0.5","W = 1","W = 1.5"),col=rainbow(N),lwd=2)
mtext(paste(np," particles",sep=""),side=3)
plot(acomp(out$mod.post[,,2]),labels=c("W = 0.5","W = 1","W = 1.5"),col=rainbow(N),lwd=2)

# Comparing W = 0.5 versus W = 1.0
mod.post2 = array(NA,dim=c(N,2,2))
mod.priors = c(.5, .5)
for(i in 1:N)
{
  mod.post2[i,,1] = out$post.lprob[i,1:2,1] + log(mod.priors)  - log(sum(exp(out$post.lprob[i,1:2,1])*mod.priors))
  mod.post2[i,,2] = out$post.lprob[i,1:2,2] + log(mod.priors)  - log(sum(exp(out$post.lprob[i,1:2,2])*mod.priors))
}
mod.post2 = exp(mod.post2)

windows(width=10,height=10)
par(mfcol=c(2,2))
r = sample(1:N, 40, rep=FALSE)
plot(rep(1,length(r)),mod.post2[r,2,1],main=paste(np," particles",sep=""),col=rainbow(N),xaxt = 'n',yaxt='n',lwd=2,xlim=c(0.9,2.1),ylim=c(0,1),xlab="",ylab=expression(paste("P(",tilde(W)," = 1) vs P(",tilde(W)," = 0.5)",sep="")))
axis(side=1, at=c(1,2), labels=c("Direct Calculation","PF approx."))
axis(side=2, at=c(0,0.5,1))
points(rep(2,length(r)),mod.post2[r,2,2],col=rainbow(length(r)),lwd=2)
segments(rep(1,length(r)), mod.post2[r,2,1], rep(2,length(r)), mod.post2[r,2,2])

# Comparing W = 1.0 versus W = 1.5
mod.post2 = array(NA,dim=c(N,2,2))
mod.priors = c(.5, .5)
for(i in 1:N)
{
  mod.post2[i,,1] = out$post.lprob[i,2:3,1] + log(mod.priors)  - log(sum(exp(out$post.lprob[i,2:3,1])*mod.priors))
  mod.post2[i,,2] = out$post.lprob[i,2:3,2] + log(mod.priors)  - log(sum(exp(out$post.lprob[i,2:3,2])*mod.priors))
}
mod.post2 = exp(mod.post2)

plot(rep(1,length(r)),mod.post2[r,1,1],col=rainbow(N),xaxt = 'n',yaxt = 'n',lwd=2,xlim=c(0.9,2.1),ylim=c(0,1),xlab="",ylab=expression(paste("P(",tilde(W)," = 1) vs P(",tilde(W)," = 1.5)",sep="")))
axis(side=1, at=c(1,2), labels=c("Direct Calculation","PF approx."))
axis(side=2, at=c(0,0.5,1))
points(rep(2,length(r)),mod.post2[r,1,2],col=rainbow(length(r)),lwd=2)
segments(rep(1,length(r)), mod.post2[r,1,1], rep(2,length(r)), mod.post2[r,1,2])

# 200 particles
np = 200
load(paste("../data/10-9-13/mc_pf_test-",nt,"-",np,"-",N,".rdata",sep=""))
# Comparing W = 0.5 versus W = 1.0
mod.post2 = array(NA,dim=c(N,2,2))
mod.priors = c(.5, .5)
for(i in 1:N)
{
  mod.post2[i,,1] = out$post.lprob[i,1:2,1] + log(mod.priors)  - log(sum(exp(out$post.lprob[i,1:2,1])*mod.priors))
  mod.post2[i,,2] = out$post.lprob[i,1:2,2] + log(mod.priors)  - log(sum(exp(out$post.lprob[i,1:2,2])*mod.priors))
}
mod.post2 = exp(mod.post2)

plot(rep(1,length(r)),mod.post2[r,2,1],main=paste(np," particles",sep=""),col=rainbow(N),xaxt = 'n',yaxt='n',lwd=2,xlim=c(0.9,2.1),ylim=c(0,1),xlab="",ylab="")
axis(side=1, at=c(1,2), labels=c("Direct Calculation","PF approx."))
axis(side=2, at=c(0,0.5,1))
points(rep(2,length(r)),mod.post2[r,2,2],col=rainbow(length(r)),lwd=2)
segments(rep(1,length(r)), mod.post2[r,2,1], rep(2,length(r)), mod.post2[r,2,2])

# Comparing W = 1.0 versus W = 1.5
mod.post2 = array(NA,dim=c(N,2,2))
mod.priors = c(.5, .5)
for(i in 1:N)
{
  mod.post2[i,,1] = out$post.lprob[i,2:3,1] + log(mod.priors)  - log(sum(exp(out$post.lprob[i,2:3,1])*mod.priors))
  mod.post2[i,,2] = out$post.lprob[i,2:3,2] + log(mod.priors)  - log(sum(exp(out$post.lprob[i,2:3,2])*mod.priors))
}
mod.post2 = exp(mod.post2)

plot(rep(1,length(r)),mod.post2[r,1,1],col=rainbow(N),xaxt = 'n',yaxt = 'n',lwd=2,xlim=c(0.9,2.1),ylim=c(0,1),xlab="",ylab="")
axis(side=1, at=c(1,2), labels=c("Direct Calculation","PF approx."))
axis(side=2, at=c(0,0.5,1))
points(rep(2,length(r)),mod.post2[r,1,2],col=rainbow(length(r)),lwd=2)
segments(rep(1,length(r)), mod.post2[r,1,1], rep(2,length(r)), mod.post2[r,1,2])

# Plot true log marginal likelihoods side by side with particle filter approximations for each model
windows(width=8, height=5)
par(mfrow=c(1,3), mar=c(5,4,4,0)+0.1)
r = sample(1:N, N, replace=FALSE)
ymin = min(out$post.lprob[r,,]); ymax = max(out$post.lprob[r,,])
plot(rep(1,length(r)), out$post.lprob[r,1,1], col=rainbow(length(r)), lwd=2, xlim = c(0.9,2.1), ylim = c(ymin,ymax), xaxt = 'n', xlab = "", ylab = "Log Marginal Likelihood", main=expression(paste(tilde(W)," = 0.5",sep="")))
axis(side=1, at=c(1,2), labels=c("Direct Calculation","PF approx."))
points(rep(2,length(r)),  out$post.lprob[r,1,2], col=rainbow(length(r)), lwd=2)
segments(rep(1,length(r)), out$post.lprob[r,1,1], rep(2,length(r)), out$post.lprob[r,1,2])
par(mar=c(5,2,4,2)+0.1)
plot(rep(1,length(r)), out$post.lprob[r,2,1], col=rainbow(length(r)), lwd=2, xlim = c(0.9,2.1), ylim = c(ymin,ymax), xaxt = 'n', xlab = "", ylab = "", main=expression(paste(tilde(W)," = 1",sep="")))
axis(side=1, at=c(1,2), labels=c("Direct Calculation","PF approx."))
points(rep(2,length(r)),  out$post.lprob[r,2,2], col=rainbow(length(r)), lwd=2)
segments(rep(1,length(r)), out$post.lprob[r,2,1], rep(2,length(r)), out$post.lprob[r,2,2])
par(mar=c(5,2,4,2)+0.1)
plot(rep(1,length(r)), out$post.lprob[r,3,1], col=rainbow(length(r)), lwd=2, xlim = c(0.9,2.1), ylim = c(ymin,ymax), xaxt = 'n', xlab = "", ylab = "", main=expression(paste(tilde(W)," = 1.5",sep="")))
axis(side=1, at=c(1,2), labels=c("Direct Calculation","PF approx."))
points(rep(2,length(r)),  out$post.lprob[r,3,2], col=rainbow(length(r)), lwd=2)
segments(rep(1,length(r)), out$post.lprob[r,3,1], rep(2,length(r)), out$post.lprob[r,3,2])

## Contingency tables
nt = 150
np = 100
N = 100
load(paste("../data/10-9-13/mc_pf_test-",nt,"-",np,"-",N,".rdata",sep=""))

# Consider two models: W = 1.0 and W = 1.5
twoMod.lprob = out$post.lprob[,2:3,]
topMod_truth = apply(twoMod.lprob[,,1], 1, function(x) which(x == max(x)))
topMod_pf = apply(twoMod.lprob[,,2], 1, function(x) which(x == max(x)))
table(topMod_truth, topMod_pf)

# Consider two models: W = 0.5 and W = 1.0
twoMod.lprob = out$post.lprob[,1:2,]
topMod_truth = apply(twoMod.lprob[,,1], 1, function(x) which(x == max(x)))
topMod_pf = apply(twoMod.lprob[,,2], 1, function(x) which(x == max(x)))
table(topMod_truth, topMod_pf)

# Consider three models: W = 0.5, 1.0, and 1.5
topMod_truth = apply(out$post.lprob[,,1], 1, function(x) which(x == max(x)))
topMod_pf = apply(out$post.lprob[,,2], 1, function(x) which(x == max(x)))
table(topMod_truth, topMod_pf)

nt = 150
np = 200
N = 100
load(paste("../data/10-9-13/mc_pf_test-",nt,"-",np,"-",N,".rdata",sep=""))

# Consider two models: W = 1.0 and W = 1.5
twoMod.lprob = out$post.lprob[,2:3,]
topMod_truth = apply(twoMod.lprob[,,1], 1, function(x) which(x == max(x)))
topMod_pf = apply(twoMod.lprob[,,2], 1, function(x) which(x == max(x)))
table(topMod_truth, topMod_pf)

# Consider two models: W = 0.5 and W = 1.0
twoMod.lprob = out$post.lprob[,1:2,]
topMod_truth = apply(twoMod.lprob[,,1], 1, function(x) which(x == max(x)))
topMod_pf = apply(twoMod.lprob[,,2], 1, function(x) which(x == max(x)))
table(topMod_truth, topMod_pf)

# Consider three models: W = 0.5, 1.0, and 1.5
topMod_truth = apply(out$post.lprob[,,1], 1, function(x) which(x == max(x)))
topMod_pf = apply(out$post.lprob[,,2], 1, function(x) which(x == max(x)))
table(topMod_truth, topMod_pf)