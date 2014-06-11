source("dlm_ar_functions.r")
source("dlm_mcmc_functions.r")
source("pl_fmri_functions.r")
source("pl.r")
source("pf_functions.r")

mod = "M101s"

# Load data
load(paste("../data/dlm_ar_sim-20-",mod,"-2.rdata",sep=""))
n = 19
n.sim = 3
mysim = get(paste(mod,"_dat",sep=""))[[1]][[n]][[n.sim]]

# Set data and design matrix
y = mysim$y
U = t(mysim$true.params$U[1,,])
nt = dim(y)[2]
truth.theta = c(mysim$true.params$beta, mysim$true.params$G[,1], mysim$true.params$W[1,1], mysim$true.params$V[1,1])
d = length(mysim$true.params$beta)
p = dim(mysim$true.params$G)[1]
same = length(strsplit(mod,"")[[1]]) == 5

# Calculate MLE fit model
require(dlm)
dyn = as.numeric(strsplit(mod,"")[[1]][2]) == 1
fit = dlmMLE(y, c(truth.theta[(d+1):(d+p)],log(truth.theta[(d+p+1):(d+p+2)])), function(par) build.ar1(par, U, dyn, same))
s=dlmSmooth(dlmFilter(y, build.ar1(fit$par, U, dyn=T)))
mle = c(s$s[nt+1,1:d],fit$par[1:p],exp(fit$par[c(p+1,p+2)]))

# Run MCMC with vague prior
psi = list(U = mysim$true.params$U, V = matrix(1), F = mysim$true.params$F)
prior = list(m0 = rep(0,p), C0 = 1e6*diag(p), b0 = rep(0,d), B0 = 1e6*diag(d), am0 = 1e-6, bm0 = 1e-6, phi0 = list(rep(0,p)), Phi0 = list(1e6*diag(p)), as0 = 1e-6, bs0 = 1e-6)
mcmc.out = dlm.ar.mcmc(y, psi, prior, mcmc.details=list(n.sims = 11000, n.burn = 1000, n.thin = 1))
x.quant = apply(mcmc.out$x[,1,], 2, function(x) quantile(x, probs=c(.025,.975)))

# Construct histograms to see if agree with MLE
windows(width=15, height=10)
par(mfrow=c(2,3), mar=c(5,7,4,2)+0.1)
burn=10
if(p == 1)
{
  ymin = min(x.quant[1,-(1:burn)], s$s[,dim(s$s)[2]])
  ymax = max(x.quant[2,-(1:burn)], s$s[,dim(s$s)[2]])
  plot(0:nt, s$s[,3], type="l", ylim=c(ymin,ymax), col=2, lwd=2, xlab=expression(t), ylab=expression(x[t]), cex.lab=1.5)
  lines(0:nt, x.quant[1,])
  lines(0:nt, x.quant[2,])
  legend("bottomright",c("MCMC","MLE"), lty=c(1,1), lwd=c(1,2), col=c(1,2))
}
hist(mcmc.out$beta[,1], xlab=expression(beta[0]), main=eval(bquote(expression(M[.(paste(strsplit(mod,"")[[1]][2:4],sep="",collapse=""))]))), cex.main=1.5)
abline(v=mle[1], col=2, lwd=2)
blabels = rep(NA, d+p+2)
truth.theta = round(truth.theta,2)
for(j in 1:d) blabels[j] = eval(bquote(expression(paste(beta[.(j-1)]," = ",.(truth.theta[j]),sep = ""))))
for(j in (d+1):(d+p)) blabels[j] = eval(bquote(expression(paste(phi[.(j-d)]," = ",.(truth.theta[j]),sep=""))))
blabels[d+p+1] = eval(bquote(expression(paste(sigma[s]^2," = ",.(truth.theta[d+p+1]),sep=""))))
blabels[d+p+2] = eval(bquote(expression(paste(sigma[m]^2," = ",.(truth.theta[d+p+2]),sep=""))))
bmin = min(hist(mcmc.out$beta[,1], plot=F)$breaks)
dist = max(hist(mcmc.out$beta[,1], plot=F)$breaks) - bmin
if(p > 1) legend("topright","MLE", lty=1, lwd=2, col=2)
for(j in 1:length(truth.theta)) mtext(blabels[j],side=3,at=bmin + ((j-1)/(d+p+2))*(dist),cex=0.8)
for(j in 2:d)
{
  hist(mcmc.out$beta[,j], xlab=eval(bquote(expression(beta[.(j-1)]))), main="")
  abline(v=mle[j], col=2, lwd=2)
}
for(j in (d+1):(d+p))
{
  hist(mcmc.out$phi[[1]][,j-d], xlab=eval(bquote(expression(phi[.(j-d)]))), main="")
  abline(v=mle[j], col=2, lwd=2)
}
hist(mcmc.out$sigma2s[,1], xlab=expression(sigma[s]^2), main="")
abline(v=mle[d+p+1], col=2, lwd=2)
hist(mcmc.out$sigma2m, xlab=expression(sigma[m]^2), main="")
abline(v=mle[d+p+2], col=2, lwd=2)

# Run MCMC with informative prior
prior = list(m0 = rep(0,p), C0 = 10*diag(p), b0 = truth.theta[1:d], B0 = 10*diag(2), am0 = 1, bm0 = 1, phi0 = list(truth.theta[(d+1):(d+p)]), Phi0 = list(.5*diag(p)), as0 = 1, bs0 = 1)
mcmc.out2 = dlm.ar.mcmc(y, psi, prior, mcmc.details=list(n.sims = 1100, n.burn = 100, n.thin = 1))
beta.quant = apply(mcmc.out2$beta, 2, function(x) quantile(x, probs=c(.025,.5,.975)))
phi.quant = apply(mcmc.out2$phi[[1]], 2, function(x) quantile(x, probs=c(.025,.5,.975)))
sigma2s.quant = apply(mcmc.out2$sigma2s, 2, function(x) quantile(x, probs=c(.025,.5,.975)))
sigma2m.quant = quantile(mcmc.out2$sigma2m, probs=c(.025,.5,.975))
mcmc.quant2 = cbind(beta.quant,phi.quant,sigma2s.quant,sigma2m.quant)
x.quant2 = apply(mcmc.out2$x, 3, function(x) quantile(x, probs=c(.025,.5,.975)))
x.cent2 = apply(mcmc.out2$x, 3, function(x) x - mean(x))
x.cent.quant2 = apply(x.cent2, 2, function(x) quantile(x, probs=c(.025,.5,.975)))

# Run pl
FF = mysim$true.params$F[1,1,]
dlpred = function(y, x, suff.x, theta) dlpred.ar(y, x, suff.x, theta, U, FF)
revo = function(y, x, suff.x, theta) revo.ar(y, x, suff.x, theta, U, FF)
smap.theta = function(suff.theta, y, x.new, x.curr) smap.theta.ar(suff.theta, y, x.new, x.curr, U, FF)
smap.state = function(suff.x, y, theta) smap.state.ar(suff.x, y, theta, U, FF)
rprior <- function(j) rprior.ar(prior)
rmove <- function(suff.theta) rmove.ar(suff.theta, prior)
out=list()
np = c(100,500,1000,5000,10000)
for(i in 1:length(np)) out[[i]] <- pl(y, dlpred, revo, rprior, rmove, smap.theta, smap.state, n=np[i], progress=TRUE, method="stratified", nonuniformity="ess", threshold=0.8, log=F)
save(out, file = paste("../data/fmri_ar_",mod,"-test.rdata",sep=""))

# Plot sequential credible intervals
windows(width=15, height=10)
par(mfrow=c(2,3), mar=c(5,7,4,2)+0.1)
ymin = rep(Inf, length(mle)+(p == 1)); ymax = rep(-Inf, length(mle)+(p == 1))
burn = 10
out.cent.quant = theta.quant = list()
for(i in 1:length(np))
{
  if(p == 1)
  {
    out.cent = apply(out[[i]]$state[1,,], 2, function(x) x - mean(x))
    out.cent.quant[[i]] = pf.quantile(array(out.cent, c(1,dim(out.cent)[1],dim(out.cent)[2])), out[[i]]$weight, function(x, param=1) x, c(.025,.5,.975))
    ymin[1] = min(ymin[1], x.cent.quant2[1,], out.cent.quant[[i]][,1,1], 0)
    ymax[1] = max(ymax[1], x.cent.quant2[3,], out.cent.quant[[i]][,1,3], 0)
  }
  theta.quant[[i]] = pf.quantile(out[[i]]$theta, out[[i]]$weight, function(x, param=1) x, c(.025,.5,.975))
  for(j in (1+(p==1)):(length(mle)+(p == 1)))
  {
    ymin[j] = min(ymin[j], theta.quant[[i]][-(1:burn),j-(p==1),1])
    ymax[j] = max(ymax[j], theta.quant[[i]][-(1:burn),j-(p==1),3])
  }
} 
if(p == 1)
{
  plot(0:nt, rep(0,nt+1), type="l", ylim=c(ymin[1],ymax[1]), col="gray70", lwd=2, xlab=expression(t), ylab=expression(x[t]), cex.lab=1.5)
  lines(0:nt, x.cent.quant2[1,])
  lines(0:nt, x.cent.quant2[2,])
  lines(0:nt, x.cent.quant2[3,])
  for(i in 1:length(np))
  {
    lines(0:nt, out.cent.quant[[i]][,1,1], col=i+1)
    lines(0:nt, out.cent.quant[[i]][,1,2], col=i+1)
    lines(0:nt, out.cent.quant[[i]][,1,3], col=i+1)
  }
  legend("topright",c(paste(np," particles",sep=""),"mcmc-smoothed","truth"), lty=c(rep(1,length(np)),1,1), col=c(2:(length(np)+1),1,"gray70"), lwd=c(rep(1,length(np)),1,2),cex=0.9)
}
expr = rep(NA, d+p+2)
for(i in 1:d) expr[i] = eval(bquote(expression(beta[.(i-1)])))
for(i in (d+1):(d+p)) expr[i] = eval(bquote(expression(phi[.(i-d)])))
expr[d+p+1] = expression(sigma[s]^2)
expr[d+p+2] = expression(sigma[m]^2)
for(i in 1:length(truth.theta))
{
  plot(0:nt, rep(truth.theta[i],nt+1), type="l", cex.lab=1.5,lwd=2, xlab=expression(t), ylab= expr[i], col="gray70", ylim = c(min(ymin[i+1], mcmc.quant2[1,i]), max(ymax[i+1], mcmc.quant2[2,i])))
  if(i == 1)
  {
    blabels = rep(NA, 5)
    for(j in 1:2) blabels[j] = eval(bquote(expression(paste(beta[.(j-1)]," = ",.(truth.theta[j]),sep = ""))))
    blabels[3] = eval(bquote(expression(paste(phi," = ",.(truth.theta[3]),sep=""))))
    blabels[4] = eval(bquote(expression(paste(sigma[s]^2," = ",.(truth.theta[4]),sep=""))))
    blabels[5] = eval(bquote(expression(paste(sigma[m]^2," = ",.(truth.theta[5]),sep=""))))
    for(j in 1:5) mtext(blabels[j],side=3,at=((j-1)/5)*(nt+1),cex=0.8)
    title(eval(bquote(expression(M[.(paste(strsplit(mod,"")[[1]][2:4],sep="",collapse=""))]))), cex.main=1.5)
    legend("bottomright",c(paste(np," particles",sep=""),"mcmc","truth"), lty=c(rep(1,length(np)),NA,1), col=c(2:(length(np)+1),1,"gray70"), lwd=c(rep(1,length(np)),1,2), pch=c(rep(NA,length(np)),"x",NA),cex=0.9)
  }
  for(j in 1:length(np))
  {
    lines(0:nt, theta.quant[[j]][,i,1],col=j+1)
    lines(0:nt, theta.quant[[j]][,i,2],col=j+1)
    lines(0:nt, theta.quant[[j]][,i,3],col=j+1)
    points(rep(nt,3), mcmc.quant2[,i], pch='x')
  }
}


