source("pf_functions.r")
source("pf_functions.r")
source("pl.r")
source("dlm_cv_functions.r")

# Set data and graphics path
dpath = "../data/"
gpath = "../graphs/"

# Load simulated data
load("../data/rw_sim.rdata")
n.sim = 1
y = mysims[[n.sim]]$y[1,]
x = mysims[[n.sim]]$x[1,]
nt = length(y)

# Functions to run rao-blackwellized particle filter
dlpred.rb = function(y,x,suff.x,theta) dlpred.ll(y,x,suff.x,theta,rb=T)
revo.rb = function(y,x,suff.x,theta) revo.ll(y,x,suff.x,theta,rb=T)

# Run particle filters and calculate quantiles for state and parameter
alpha = 0.05

# Calculate true 95% CI of states, one-step ahead predictions, and precision
post = cv.post(y, F=1, G=1, V=1, W=1, a0=1, b0=1, m0=0, C0=1)
lk = qt(alpha/2,2*post$a)*sqrt(post$C[1,1,]*(post$b/post$a)) + post$m[1,]
uk = qt(1 - alpha/2,2*post$a)*sqrt(post$C[1,1,]*(post$b/post$a)) + post$m[1,]
lt = qt(alpha/2,2*post$a[1:nt])*sqrt(post$Q[1,1,]*(post$b[1:nt]/post$a[1:nt])) + post$f[1,]
ut = qt(1 - alpha/2,2*post$a[1:nt])*sqrt(post$Q[1,1,]*(post$b[1:nt]/post$a[1:nt])) + post$f[1,]
lp = qgamma(alpha/2,post$a,post$b)
up = qgamma(1-alpha/2,post$a,post$b)

# Plot true filtered distributions of state and precision
pdf(paste(gpath,"cv-states-precision-",n.sim,".pdf",sep=""),width=10,height=5)
par(mfrow=c(1,2), mar=c(5,6,4,2)+0.1)
ymax = max(uk,ut,x,y)
ymin = min(lk,lt,x,y)
plot(0:nt,x,ylim=c(ymin,ymax),type="l",xlab=expression(t),ylab=expression(x[t]),main="95% C.I. for filtered states",cex.lab=1.75,cex.main=1.75)
points(1:nt,y)
lines(0:nt,uk,col=2);lines(0:nt,lk,col=2)
lines(1:nt,lt,col=4);lines(1:nt,ut,col=4)
legend("topleft",c(expression(x[t]),expression(y[t]),"95% CI","95% PI"),cex=0.85,lty=c(1,NA,1,1),pch=c(NA,1,NA,NA),col =c(1,1,2,4),bg="white")
ymax = max(up,1)
ymin = min(lp,1)
plot(0:nt,up,ylim=c(ymin,ymax),type="l",col=2, xlab=expression(t),ylab=expression(1/theta),main="95% C.I. for filtered precision",cex.lab=1.75,cex.main=1.75)
lines(0:nt,lp,col=2)
abline(h=1)
legend("topright",c("95% CI","Truth"),lty=c(1,1),col=c(2,1),bg="white",cex=0.85)
dev.off()

# Plot true log marginal likelihood for increasing SNR
require(plyr)
true_lmarglik = function(lambda)
{
  post = cv.post(mysims[[n.sim]]$y, F=1, G=1, V=1, W=lambda, a0=1, b0=1, m0=0, C0=1)
  return(cv.lmarglik(mysims[[n.sim]]$y, post$f, post$Q, post$a, post$b))
}
true.lmarglik = maply(data.frame(lambda=seq(0.1,10,.1)), true_lmarglik)
windows()
par(mar=c(5,6,4,2)+0.1)
plot(seq(0.1,10,.1),true.lmarglik,type="l",xlab=expression(lambda),ylab=expression(log(p(y[1:T]))),main=expression(paste(log(p(y[1:T]))," vs ",lambda,sep="")),cex.lab=1.75,cex.main=1.75)
abline(v=c(0.5,1,2),lty=c(2,1,2))
mtext(c(0.5,1,2),side=1,at=c(0.5,1,2))

# Test particle learning
np = 5000

# Regular pl
out.pl = pl(y, dlpred.ll, revo.ll, function(j) rprior.ll(), rmove.ll, smap.theta.ll, smap.state.ll, np, progress=TRUE, method="stratified", nonuniformity="none", log=F) # always resample
state.quant.pl = pf.quantile(out.pl$state, out.pl$weight, function(x, param=1) x, c(alpha/2,1-alpha/2)) 
# theta.quant.pl = pf.quantile(1/out.pl$theta, out.pl$weight, function(x, param=1) x, c(alpha/2,1-alpha/2))
theta.quant.pl = pf.mix.quantile(list(out.pl$suff.theta), out.pl$weight, list(function(x, o, w) sum(w*pgamma(x,o[1,],o[2,])))) # Calculate quantiles from mixture distribution

# Rao-Blackwellized
out.rb = pl(y, dlpred.rb, revo.rb, function(j) rprior.ll(), rmove.ll, smap.theta.ll, smap.state.ll, np, progress=TRUE, method="stratified", nonuniformity="none", log=F) # always resample
# state.quant.rb = pf.quantile(out.rb$state, out.rb$weight, function(x, param=1) x, c(alpha/2,1-alpha/2))
# theta.quant.rb = pf.quantile(1/out.rb$theta, out.rb$weight, function(x, param=1) x, c(alpha/2,1-alpha/2))
state.quant.rb = pf.mix.quantile(list(out.rb$suff.state), out.rb$weight, list(function(x, o, w) sum(w*pnorm(x,o[1,],sqrt(o[2,])))))
theta.quant.rb = pf.mix.quantile(list(out.rb$suff.theta), out.rb$weight, list(function(x, o, w) sum(w*pgamma(x,o[1,],o[2,]))))

# Plot 95% CI for states and precision
windows(width=10,height=5)
par(mfrow=c(1,2))
nt = dim(y)[2]
gmin = min(mysims[[n.sim]]$x[1,], lk, state.quant.rb[,1,1], state.quant.pl[,1,1])
gmax = max(mysims[[n.sim]]$x[1,], uk, state.quant.rb[,1,2], state.quant.pl[,1,2])
plot(0:nt,mysims[[n.sim]]$x[1,],ylim=c(gmin,gmax),type="l",xlab=expression(t),ylab="Position",col='gray',lwd=3)
mtext(paste(np," particles",sep=""),side=3)
lines(0:nt,lk,lwd=3)
lines(0:nt,uk,lwd=3)
lines(0:nt, state.quant.rb[,1,1],col=4)
lines(0:nt, state.quant.rb[,1,2],col=4)
lines(0:nt, state.quant.pl[,1,1],col=2)
lines(0:nt, state.quant.pl[,1,2],col=2)
legend("topleft",legend=c('pl','pl-rb', "True Post", "True Sim"),lty=c(1,1,1,1),col=c(2,4,1,'gray'),lwd=c(1,1,3,3))
title("95% CI for filtered states")

gmin = min(1, lp, theta.quant.rb[,1,1], theta.quant.pl[,1,1])
gmax = max(1, up, theta.quant.rb[,1,2], theta.quant.pl[,1,2])
plot(0:nt,lp,type="l",lwd=3,ylim=c(gmin,gmax),xlab=expression(t),ylab=expression(1/theta))
mtext(paste(np," particles",sep=""),side=3)
lines(0:nt,lp,lwd=3)
lines(0:nt,up,lwd=3)
lines(0:nt, theta.quant.rb[,1,1],col=4)
lines(0:nt, theta.quant.rb[,1,2],col=4)
lines(0:nt, theta.quant.pl[,1,1],col=2)
lines(0:nt, theta.quant.pl[,1,2],col=2)
abline(h=1,col='gray',lwd=3)
legend("topright",legend=c('pl','pl-rb', "True Post", "True Sim"),lty=c(1,1,1,1),col=c(2,4,1,'gray'),lwd=c(1,1,3,3))
title("95% CI for filtered precision")