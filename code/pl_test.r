source("pf_functions.r")
source("pl_ll_functions.r")
source("pl.r")
source("dlm_cv_functions.r")

# Load simulated data
load("../data/rw_sim.rdata")
n.sim = 1
y = mysims[[n.sim]]$y

# Functions to run rao-blackwellized particle filter
dlpred.rb = function(y,x,suff.x,theta) dlpred.ll(y,x,suff.x,theta,rb=T)
revo.rb = function(y,x,suff.x,theta) revo.ll(y,x,suff.x,theta,rb=T)

# Run particle filters and calculate quantiles for state and parameter
np = 5000
alpha = 0.05

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

# Calculate true 95% CI of states and precision
post = cv.post(y, F=1, G=1, V=1, W=1, a0=1, b0=1, m0=0, C0=1)
lk = qt(alpha/2,2*post$a)*sqrt(post$C[1,1,]*(post$b/post$a)) + post$m[1,]
uk = qt(1 - alpha/2,2*post$a)*sqrt(post$C[1,1,]*(post$b/post$a)) + post$m[1,]
lp = qgamma(alpha/2,post$a,post$b)
up = qgamma(1-alpha/2,post$a,post$b)

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