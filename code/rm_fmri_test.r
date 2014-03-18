source("rm_pf.r")
source("pf_functions.r")
source("pf_mc_functions.r")

# Set data and graphics path
dpath = "../data/"
gpath = "../graphs/"

# Test dlm MCMC functions
source("rm_fmri_functions.r")

require(dlm)
N = 20
n = 6
n.sim = 1
mod = "M101"
np = 100
  
# Load data
load(paste(dpath,"dlm_ar_sim-",N,"-",mod,".rdata",sep=""))
mysim = get(paste(mod,"_dat",sep=""))[[1]][[n]][[n.sim]]
y = mysim$y
  
# Set known values and get dimensions of beta and phi
U = mysim$true.params$U
F = mysim$true.params$F
V = mysim$true.params$V
psi = list(U=U,F=F,V=V)
d = dim(U)[2]
p = dim(F)[2]
nt = dim(U)[3]
  
# Calculate MLEs
fit = lm(y[1,] ~ t(U[1,,]) - 1)
phi.init = as.numeric(acf(residuals(fit),type='partial',plot=FALSE)$acf)[1:p]
sigma2m.init = summary(fit)$sigma^2 / 2
sigma2s.init = summary(fit)$sigma^2 / 2
build <- function(par) get(paste("build_",mod,sep=""))(par,V,U,U[1,2,])
fit.mle <- dlmMLE(y, c(logit(phi.init,-1,1),log(sigma2s.init),log(sigma2m.init)),build, hessian = TRUE)
par.var <- solve(fit.mle$hessian)
fit.smooth <- dlmSmooth(dlmFilter(y, build(fit.mle$par)))
svar = dlmSvd2var(fit.smooth$U.S, fit.smooth$D.S)
  
x.mle = t(fit.smooth$s[,(d+1):(d+p)])
x0.var = svar[[nt+1]][d+1,d+1]
beta.mle = fit.smooth$s[nt+1,1:d]
beta.var = svar[[nt+1]][1:d,1:d]
phi.mle = unlogit(fit.mle$par[1:p],-1,1)
if(p == 1) phi.var = (d.unlogit(fit.mle$par[1:p],-1,1)^2) * par.var[1,1] else phi.var = diag(p)
sigma2m.mle = exp(fit.mle$par[p+2])
sigma2m.var = exp(fit.mle$par[p+2]) * par.var[p+2,p+2]
sigma2s.mle = exp(fit.mle$par[p+1])
sigma2s.var = exp(fit.mle$par[p+1]) * par.var[p+1,p+1]
theta.mle = c(beta.mle,phi.mle,sigma2m.mle,sigma2s.mle)
theta.truth = c(mysim$true.params$beta,mysim$true.params$G[,1],mysim$true.params$V[1,1],mysim$true.params$W[1,1])

# Load MCMC data and calculate 95% CI from MCMC samples
out.all <- list()
for(i in 1:3)
{
  load(paste(dpath,"fmri_dlm_mcmc_test-",paste(N,n,n.sim,mod,i,11000,1000,10,1,1,1,1,1,sep="-"),".rdata",sep=""))
  out.all[[i]] = out
}
x.mcmc <- quantile(sapply(out.all, function(a) a$x[,1,nt+1]), c(0.025,.975))
beta0.mcmc <- quantile(sapply(out.all, function(x) x$beta[,1]), c(.025,.975))
beta1.mcmc <- quantile(sapply(out.all, function(x) x$beta[,2]), c(.025,.975))
phi.mcmc <- quantile(sapply(out.all, function(x) x$phi[,1]), c(.025,.975))
sigma2m.mcmc <- quantile(sapply(out.all, function(x) x$sigma2m), c(.025,.975))
sigma2s.mcmc <- quantile(sapply(out.all, function(x) x$sigma2s), c(.025,.975))
theta.mcmc <- rbind(beta0.mcmc,beta1.mcmc,phi.mcmc,sigma2m.mcmc,sigma2s.mcmc)

# Define MCMC kernel, priors, likelihood, and run resample-move particle filter
dllik <- function(y, x, theta) dllik.dlm(y, x, theta, t(U[1,,]), U[1,2,])
prior = list(b0 = beta.mle, B0 = beta.var, phi0 = phi.mle, Phi0 = phi.var, am0 = sigma2m.mle^2 / sigma2m.var, bm0 = sigma2m.mle / sigma2m.var , as0 = sigma2s.mle^2 / sigma2s.var, bs0 = sigma2s.mle / sigma2s.var, m0 = rep(0,p))
rprior <- function() rprior.dlm(prior$b0, prior$B0, prior$phi0, prior$Phi0, prior$am0, prior$bm0, prior$as0, prior$bs0, prior$m0)
rmove <- function(y, x, theta) rmove.dlm(y, x, theta, psi, prior)
out = rm_pf(y, dllik, revo.dlm, rprior, rmove, np, method="stratified", nonuniformity="ess", threshold=0.8, log=FALSE)

# Calculate approximated 95% CI of states and precision
alpha = 0.05
state.quant = pf.quantile(out$state, out$weight, function(x, param=1) x, c(alpha/2,1-alpha/2))
theta.quant = pf.quantile(out$theta, out$weight, function(x, param=1) x, c(alpha/2,1-alpha/2))
  
# Plot 95% CI and compare with MLEs, truth, MCMC 95% CI
nt = dim(out$state)[3] - 1
burn = 1
pdf(paste(gpath,"rm-fmri-test-",paste(N, n, n.sim, mod, np,sep="-"),".pdf",sep=""),width=15,height=10)
par(mfrow = c(2,3),mar=c(5,7,4,2)+0.1)
gmin = min(mysim$x, x.mle[1,], x.mcmc[1], state.quant[-(1:burn),1,1])
gmax = max(mysim$x, x.mle[1,], x.mcmc[2], state.quant[-(1:burn),1,2])
plot(0:nt,state.quant[,1,1],type="l",ylim=c(gmin,gmax),ylab=expression(x[t]),xlab=expression(t),cex.lab=2,cex.axis=1.3)
lines(0:nt,state.quant[,1,2])
lines(0:nt,x.mle[1,],col="gray")
lines(0:nt,mysim$x,col=2)
abline(h=x.mcmc, col = 4)
ylabs = expression(beta[0],beta[1],phi[1],sigma[m]^2,sigma[s]^2)
for(i in 1:5)
{
  if(i == 3) burn = 12 else burn = 1
  gmin = min(theta.truth[i], theta.mle[i], theta.mcmc[i,1], theta.quant[-(1:burn),i,1])
  gmax = max(theta.truth[i], theta.mle[i], theta.mcmc[i,2], theta.quant[-(1:burn),i,2])
  plot(0:nt,theta.quant[,i,1],type="l",ylim=c(gmin,gmax),ylab=ylabs[i],xlab="",cex.lab=2,cex.axis=1.3)
  lines(0:nt,theta.quant[,i,2])
  abline(h=theta.mle[i],col="gray")
  abline(h=theta.truth[i],col=2)
  abline(h=theta.mcmc[i,], col = 4)
  if(i == 3) legend("bottomright",c("95% CI (RM)","95% CI (MCMC)","MLE","Truth"),lty=c(1,1,1,1),col=c(1,4,"gray",2),cex=1.5)
}
dev.off()
  
# Calculate log marginal likelihood
pf.lmarglik(out)

# Save RM data
save(out, file = paste(dpath,"rm-fmri-test-",paste(N,n,n.sim,mod,np,sep="-"),".rdata",sep=""))
