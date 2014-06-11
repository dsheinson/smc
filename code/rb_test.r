# Test particle filters on local level model
source("pf_mc_functions.r")
source("dlm_cv_functions.r")
source("pf_functions.r")
source("pl_ll_functions.r")
source("pl.r")
require(smcUtils)

# Set data and graphics path
dpath = "../data/"
gpath = "../graphs/"

# Set simulation parameters
lambda=1
np=10000
n.sim=1
alpha=0.05

# Load simulated data
load(paste(dpath,"rw_sim.rdata",sep=""))
y = mysims[[n.sim]]$y
nt = dim(y)[2]

# Calculate true 95% CI of states and precision
post = cv.post(mysims[[n.sim]]$y, F=1, G=1, V=1, W=lambda, a0=1, b0=1, m0=0, C0=1)
lk = qt(alpha/2,2*post$a)*sqrt(post$C[1,1,]*(post$b/post$a)) + post$m[1,]
uk = qt(1 - alpha/2,2*post$a)*sqrt(post$C[1,1,]*(post$b/post$a)) + post$m[1,]
lp = qgamma(alpha/2,post$a,post$b)
up = qgamma(1-alpha/2,post$a,post$b)

##############################################################################

# Load prior samples and compute weights conditional and states and state sufficient statistics
x = rep(NA, np)
w1 = w2 = rep(NA, np)
suff.m = suff.v = suff.a = suff.b = rep(NA, np)
theta = rep(NA, np)
for(i in 1:np)
{
  tmp = rprior.ll()
  x[i] = tmp$x[1]
  theta[i] = tmp$theta
  suff.m[i] = tmp$suff.x[1]
  suff.v[i] = tmp$suff.x[2]
  suff.a[i] = tmp$suff.theta[1]
  suff.b[i] = tmp$suff.theta[2]
  w1[i] = dlpred.ll(y[,1], x[i], tmp$suff.x, tmp$theta)
  w2[i] = dlpred.ll(y[,1], x[i], tmp$suff.x, tmp$theta, rb=T)
}
w1 = renormalize(w1,log=T) # renormalize weights to [0,1]
w2 = renormalize(w2,log=T)
x11()
plot(density(log(w1)),type="l",xlim=c(-20,0)) # Are weights based on suff. stats flatter?
lines(density(log(w2)),col=2)

# Resample particles according each weight distribution
ind1 = resample(w1, method="stratified", nonunif="none", log=F)$indices
ind2 = resample(w2, method="stratified", nonunif="none", log=F)$indices
r1 = x[ind1] # resampled states
r2 = x[ind2]
suff.mr1 <- suff.m[ind1] # suff stats for states
suff.vr1 <- suff.v[ind1]
suff.mr2 <- suff.m[ind2] 
suff.vr2 <- suff.v[ind2]
theta1 = theta[ind1] # theta
theta2 = theta[ind2]
suff.a1 = suff.a[ind1] # suff stats for theta
suff.b1 = suff.b[ind1]
suff.a2 = suff.a[ind2]
suff.b2 = suff.b[ind2]

# Diagnostic plots
windows(width=15,height=10)
par(mfcol=c(2,3))
plot(density(log(theta)),type="l") # Are suff stats for theta accurate at time 0?
for(i in 1:1000) lines(density(log(1/rgamma(np,suff.a,suff.b))),col=2)
lines(density(log(theta)),lwd=2)
plot(density(x),type="l",xlim=c(-10,10),ylim=c(0,.5)) # Are suff. stats for state accurate?
for(i in 1:1000) lines(density(rnorm(np,suff.m,sqrt(suff.v))),col=2)
lines(density(x),lwd=2)
plot(density(log(theta1)),type="l") # Are resampled suff stats for theta accurate?
for(i in 1:1000) lines(density(log(1/rgamma(np,suff.a1,suff.b1))),col=2)
lines(density(log(theta1)),lwd=2)
lines(density(log(theta)),lty=2,lwd=2)
plot(density(r1),type="l",xlim=c(-10,10),ylim=c(0,.5)) # Are resampled suff. stats for state accurate?
for(i in 1:1000) lines(density(rnorm(np,suff.mr1,sqrt(suff.vr1))),col=2)
lines(density(r1),lwd=2)
lines(density(x),lty=2,lwd=2)
plot(density(log(theta2)),type="l") # Are resampled suff stats for theta accurate?
for(i in 1:1000) lines(density(log(1/rgamma(np,suff.a2,suff.b2))),col=2)
lines(density(log(theta2)),lwd=2)
lines(density(log(theta)),lty=2,lwd=2)
plot(density(r2),type="l",xlim=c(-10,10),ylim=c(0,.5)) # Are resampled suff. stats for state accurate?
for(i in 1:1000) lines(density(rnorm(np,suff.mr2,sqrt(suff.vr2))),col=2)
lines(density(r2),lwd=2)
lines(density(x),lty=2,lwd=2)

# Propagate particles and suff stats using resampled states opposed to resampled suff stats
x11 = x12 = x21 = x22 = xs = matrix(NA, 2, np)
suffb11 = suffb12 = suffb21 = suffb22 = rep(NA, np)
suff.s = matrix(NA, 2, np)
suff.bs = rep(NA, np)
theta.s = theta11 = theta12 = theta21 = theta22 = rep(NA, np)
for(i in 1:np)
{
  x11[,i] = revo.ll(y[,1], r1[i], c(suff.mr1[i],suff.vr1[i]), theta1[i])
  x12[,i] = revo.ll(y[,1], r1[i], c(suff.mr1[i],suff.vr1[i]), theta1[i], rb=T)
  x21[,i] = revo.ll(y[,1], r2[i], c(suff.mr2[i],suff.vr2[i]), theta2[i])
  x22[,i] = revo.ll(y[,1], r2[i], c(suff.mr2[i],suff.vr2[i]), theta2[i], rb=T)
  
  suffb11[i] = smap.theta.ll(c(suff.a1[i],suff.b1[i]), y[,1], x11[,i])[2]
  suffb12[i] = smap.theta.ll(c(suff.a1[i],suff.b1[i]), y[,1], x12[,i])[2]
  suffb21[i] = smap.theta.ll(c(suff.a2[i],suff.b2[i]), y[,1], x21[,i])[2]
  suffb22[i] = smap.theta.ll(c(suff.a2[i],suff.b2[i]), y[,1], x22[,i])[2]

  suff.s[,i] = smap.state.ll(c(suff.mr2[i],suff.vr2[i]), y[,1], theta2[i])
  xs[,i] = c(rnorm(1, suff.s[1,i], sqrt(suff.s[2,i])), r2[i])
  suff.bs[i] = smap.theta.ll(c(suff.a2[i],suff.b2[i]), y[,1], xs[,i])[2]
  
  theta11[i] = 1 / rgamma(1, 2.5, suffb11[i])
  theta12[i] = 1 / rgamma(1, 2.5, suffb12[i])
  theta21[i] = 1 / rgamma(1, 2.5, suffb21[i])
  theta22[i] = 1 / rgamma(1, 2.5, suffb22[i])
  theta.s[i] = 1 / rgamma(1, 2.5, suff.bs[i])
}
x11()
plot(density(x22[1,]),type="l",xlim=c(-10,10))
for(i in 1:1000) lines(density(rnorm(1000,suff.s[1,],sqrt(suff.s[2,]))),col=2)

# Diagnostics
windows(width=10,height=10)
par(mfrow=c(2,2))
plot(density(log(suffb11)),type="l") # distribution of suff stats for theta
lines(density(log(suffb12)),lty=2)
lines(density(log(suffb21)),col=2)
lines(density(log(suffb22)),col=2,lty=2)
lines(density(log(suff.bs)),col=4)
plot(density(x11[2,]),type="l",xlim=c(-5,5)) # distribution of current states
lines(density(x12[2,]),lty=2)
lines(density(x21[2,]),col=2)
lines(density(x22[2,]),lty=2,col=2)
lines(density(xs[2,]),col=4)
lines(density(r1),col=3)
lines(density(r2),col=3,lty=2)
abline(v=y[,1])
plot(density(x11[1,]),type="l",xlim=c(-5,5)) # distribution of propagated states
lines(density(x12[1,]),lty=2)
lines(density(x21[1,]),col=2)
lines(density(x22[1,]),lty=2,col=2)
lines(density(xs[1,]),col=4)
abline(v=y[,1])
plot(density(theta11),type="l",xlim=c(0,40))
lines(density(theta12),lty=2)
lines(density(theta21),col=2)
lines(density(theta22),lty=2,col=2)
lines(density(theta.s),col=4)
abline(v=c(lp[2],up[2]))
