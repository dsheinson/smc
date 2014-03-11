source("dlm_cv_functions.r")
source("pf_mc_functions.r")
source("pf_functions.r")

# Set data, graphics, and tables paths
dpath = "../data/"
gpath = "../graphs/"
tpath = "../tables/"

# Load simulated data
load(paste(dpath,"rw_sim.rdata",sep=""))
N = length(mysims)
np = c(500, 1000, 5000, 10000)
cols = rainbow(length(np))

## Function to plot quantiles from rm_pf, compared against true posterior
cv_pf_quantiles <- function(n.sim, nruns, W, alpha = 0.05, burn.k = 1, burn.p = 5)
{
  F = mysims[[n.sim]]$true.params$F[1,1,1]
  G = mysims[[n.sim]]$true.params$G[1,1]
  V = mysims[[n.sim]]$true.params$V[1,1]
  m0 = 0
  C0 = 1
  a0 = b0 = 1
  nt = dim(mysims[[n.sim]]$y)[2]
  
  # Calculate true posterior distribution under correct model
  post = cv.post(mysims[[n.sim]]$y, F, G, V, W, a0, b0, m0, C0)
  
  # Calculatue 95% credible intervals for filtered states, precision, and one-step ahead predictions
  lk = qt(alpha/2,2*post$a)*sqrt(post$C[1,1,]*(post$b/post$a)) + post$m[1,]
  uk = qt(1 - alpha/2,2*post$a)*sqrt(post$C[1,1,]*(post$b/post$a)) + post$m[1,]
  lp = qgamma(alpha/2,post$a,post$b)
  up = qgamma(1-alpha/2,post$a,post$b)
  
  # Load quantiles from rm pf to find min and max values
  gmin.k = gmin.p = Inf
  gmax.k = gmax.p = -Inf
  for(j in 1:length(np))
  {
    for(i in 1:nruns)
    {
      load(paste(dpath,"cv_pf-",paste(n.sim, i, W, np[j], alpha, sep="-"),".rdata",sep=""))
      gmin.k = min(gmin.k, lk[-(1:burn.k)], mysims[[n.sim]]$x[1,], pf.out$state.quant[-(1:burn.k),1,1])
      gmax.k = max(gmax.k, uk[-(1:burn.k)], mysims[[n.sim]]$x[1,], pf.out$state.quant[-(1:burn.k),1,2])
      gmin.p = min(gmin.p, lp[-(1:burn.p)], 1, pf.out$theta.quant[-(1:burn.p),1,1])
      gmax.p = max(gmax.p, up[-(1:burn.p)], 1, pf.out$theta.quant[-(1:burn.p),1,2])
    }
  }
  
  # Plot true data, states, 95% CI of filtered states, and 95% one-step ahead PI
  pdf(paste(gpath,"cv-pf-states-",n.sim,"-",W,"-",alpha,".pdf",sep=""))
  for(j in 1:length(np))
  {
    for(i in 1:length(nruns))
    {
      load(paste(dpath,"cv_pf-",paste(n.sim, i, W, np[j], alpha, sep="-"),".rdata",sep=""))
      if(j == 1 & i == 1) plot(0:nt,mysims[[n.sim]]$x[1,],ylim=c(gmin.k,gmax.k),type="l",xlab=expression(t),ylab=expression(x[t]),cex.lab=1.5)
      lines(0:nt, pf.out$state.quant[,1,1],col=cols[j])
      lines(0:nt, pf.out$state.quant[,1,2],col=cols[j])
    }
  }
  lines(0:nt,lk,col=2,lwd=4)
  lines(0:nt,uk,col=2,lwd=4)
  title(expression(paste("Sequential 95% CI for ",sep=""),x[t]))
  mtext(substitute(paste(tilde(W)," = ",aa),list(aa = W)),side=3)
  legend("bottomright", c(paste(np," particles",sep=""),"True post", "True sim"), lwd=c(rep(1,length(np)),4,1), lty=rep(1,length(np)+2), col = c(cols,2,1))
  dev.off()
  
  # Plot 95% CI of filtered precision
  pdf(paste(gpath,"cv-pf-precision-",n.sim,"-",W,"-",alpha,".pdf",sep=""))
  for(j in 1:length(np))
  {
    for(i in 1:nruns)
    {
      load(paste(dpath,"cv_pf-",paste(n.sim, i, W, np[j], alpha, sep="-"),".rdata",sep=""))
      if(j == 1 & i == 1) plot(0:nt,rep(1,nt+1),ylim=c(gmin.p,gmax.p),type="l",xlab=expression(t),ylab=expression(phi),cex.lab=1.5)
      lines(0:nt, pf.out$theta.quant[,1,1],col=cols[j])
      lines(0:nt, pf.out$theta.quant[,1,2],col=cols[j])
    }
  }
  lines(0:nt,lp,col=2,lwd=4)
  lines(0:nt,up,col=2,lwd=4)
  title(expression(paste("Sequential 95% CI for ",sep=""),phi))
  mtext(substitute(paste(tilde(W)," = ",aa),list(aa = W)),side=3)
  legend("topright", c(paste(np," particles",sep=""),"True post", "True sim"), lwd=c(rep(1,length(np)),4,1), lty=rep(1,length(np)+2), col = c(cols,2,1))
  dev.off()
}

require(plyr)
mydata = expand.grid(n.sim = 1, nruns = 20, W = c(0.1,0.5,1,2,3))
m_ply(mydata, cv_pf_quantiles)

## Plot kernel density estimates of log-likelihood under each model, ternary compositional plot, binary plots
cv_pf_loglik <- function(n.sim, nruns, alpha = 0.05)
{
  F = mysims[[n.sim]]$true.params$F[1,1,1]
  G = mysims[[n.sim]]$true.params$G[1,1]
  V = mysims[[n.sim]]$true.params$V[1,1]
  m0 = 0
  C0 = 1
  a0 = b0 = 1
  nt = dim(mysims[[n.sim]]$y)[2]

  # Calculate true log marginal likelihoods under each model
  true_lmarglik = function(W)
  {
    post = cv.post(mysims[[n.sim]]$y, F, G, V, W, 1, 1, 0, 1)
    return(cv.lmarglik(mysims[[n.sim]]$y, post$f, post$Q, post$a, post$b))
  }
  true.lmarglik = maply(data.frame(W=W), true_lmarglik)

  # Load approximate log marginal likelihoods for particle filter runs
  rm.lmarglik = array(NA, c(nruns, length(np), length(W)))
  for(j in 1:length(np))
  {
    for(k in 1:length(W))
    {
      for(i in 1:nruns)
      {
        load(paste(dpath,"cv_pf-",paste(n.sim, i, W[k], np[j], alpha,sep="-"),".rdata",sep=""))
        rm.lmarglik[i,j,k] = pf.out$lmarglik
      }
    }
  }

  # Plot kernel density estimates of of log marginal likelihoods for pfs compared with truth
  max.densities = aaply(rm.lmarglik, 2:3, function(a) max(density(a)$y,na.rm=TRUE))
  min.densities = aaply(rm.lmarglik, 2:3, function(a) min(density(a)$y,na.rm=TRUE))
  max.breaks = aaply(rm.lmarglik, 2:3, function(a) max(density(a)$x,na.rm=TRUE))
  min.breaks = aaply(rm.lmarglik, 2:3, function(a) min(density(a)$x,na.rm=TRUE))
  ylabs = c("Density","","")
  xlabs = c("Log Marginal Likelihood","","")
  cols = rainbow(length(np))
  mains = c(substitute(paste(tilde(W)," = ",aa,sep=""),list(aa=W[1])),substitute(paste(tilde(W)," = ",aa,sep=""),list(aa=W[2])),substitute(paste(tilde(W)," = ",aa,sep=""),list(aa=W[3])))
  pdf(file=paste(gpath,"cv_pf_loglik-",n.sim,".pdf",sep=""),width=5*length(W),height=5)
  par(mfrow=c(1,length(W)),mar=c(5,6,4,2)+0.1)
  for(i in 1:length(W))
  {
    plot(density(rm.lmarglik[,1,i]),lwd=2,col=cols[1],main=mains[i],xlab=xlabs[i],ylab=ylabs[i],xlim=c(min(min.breaks),max(max.breaks)),ylim=c(min(min.densities[,i]),max(max.densities[,i])),cex.axis=1.5,cex.lab=1.75,cex.main=2)
    if(length(np > 1)) for(j in 2:length(np)) lines(density(rm.lmarglik[,j,i]),lwd=2,col=cols[j])
    abline(v=true.lmarglik[i],lwd=1,col=1)
    for(j in which(!(1:length(W) %in% i))) abline(v=true.lmarglik[j],lwd=1,col=1,lty=2)
    mtext(round(true.lmarglik[i],2),side=1,at=true.lmarglik[i],cex=0.75)
    if(i == 1) legend("topleft",c(paste(np," particles",sep=""),"Truth","Truth (other models)"),lty=c(rep(1,length(np)),1,2),lwd=c(rep(2,length(np)),1,1),col=c(cols,1,1),cex=1.5)
  }
  dev.off()

  # Compute true posterior model probabilities
  true.postModProbs = postModProbs(true.lmarglik, rep(1/length(W), length(W)))
  
  # Compute rm_pf approximate posterior model probabilities
  rm.postModProbs = aaply(rm.lmarglik, 1:2, function(x) postModProbs(x, rep(1/length(W), length(W))))
  
  # How often does the pf pick the "right" model?
  rm.max = apply(rm.postModProbs, 1:2, function(x) which(x == max(x)))
  prob.matchTruth = apply(rm.max, 2, function(x) sum(x == which(true.postModProbs == max(true.postModProbs)))/nruns)
  
  # Ternary diagrams of posterior model probabilities
  require(compositions)
  wlabels = rep(NA,length(W))
  for(i in 1:length(W)) wlabels[i] = as.expression(bquote(paste(tilde(W)," = ",sep="") ~ .(W[i])))
  pdf(file=paste(gpath,"cv_pf_ternary-",n.sim,".pdf",sep=""),width=10,height=10)
  par(mfrow=c(2,2))
  for(i in 1:length(np))
  {
    plot(acomp(rm.postModProbs[,i,]),labels=wlabels,lwd=2)
    plot(acomp(true.postModProbs),add=TRUE,col=2,lwd=2)
    mtext(paste(np[i]," particles",sep=""),side=3,cex=2)
    if(i==1) legend("topleft",legend="True Posterior",pch=1,col=2,pt.lwd=2)
  }
  dev.off()
}

m_ply(data.frame(n.sim = 1, nruns = 20),cv_pf_loglik)

## Plots analyzing rm pf runs between simulations
# Calculate true log marginal likelihoods under each model

true_lmarglik = function(n.sim, W)
{
  F = mysims[[n.sim]]$true.params$F[1,1,1]
  G = mysims[[n.sim]]$true.params$G[1,1]
  V = mysims[[n.sim]]$true.params$V[1,1]
  m0 = 0
  C0 = 1
  a0 = b0 = 1
  post = cv.post(mysims[[n.sim]]$y, F, G, V, W, a0, b0, m0, C0)
  return(cv.lmarglik(mysims[[n.sim]]$y, post$f, post$Q, post$a, post$b))
}
true.lmarglik = maply(expand.grid(n.sim=n.sim, W=W), true_lmarglik)

# Load approximate log marginal likelihoods for particle filter runs
alpha = 0.05
rm.lmarglik = array(NA, c(length(nsims),length(np),length(W)))
for(j in 1:length(np))
{
  for(k in 1:length(W))
  {
    for(i in N)
    {
      load(paste(dpath,"cv_pf-",paste(i, 1, W[k], np[k], alpha, sep="-"),".rdata",sep=""))
      rm.lmarglik[i,j,k] = pf.out$lmarglik
    }
  }
}

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
  print(xtable(tab,caption=paste("Contingency table of models with highest posterior model probability among three models calculated analytically and approximated using the resample-move particle filter with ",np[i]," particles.",sep=""),label=paste("tab:cont-",np[i],sep="")),file=paste(tpath,"cv_pf_contTable-",np[i],".txt",sep=""),type="latex") 
}

# Plot ternary diagrams for posterior model probabilities
require(compositions)
N=length(nsims)
grad = order(true.postModProbs[,2],decreasing=TRUE)
wlabels = rep(NA,length(W))
for(i in 1:length(W)) wlabels[i] = as.expression(bquote(paste(tilde(W)," = ",sep="") ~ .(W[i])))
pdf(file=paste(gpath,"cv_pf_ternary-bw.pdf",sep=""),width=10,height=10)
par(mfrow=c(2,2))
plot(acomp(true.postModProbs[grad,]),labels=wlabels,col=gray.colors(N),lwd=2)
mtext("True Posterior",side=3,cex=2)
for(i in 2:length(np))
{
  plot(acomp(rm.postModProbs[grad,i,]),labels=wlabels,col=gray.colors(N),lwd=2)
  mtext(paste(np[i]," particles",sep=""),side=3,cex=2)
}
dev.off()

# Plot binary diagrams for posterior model probabilities
pdf(file=paste(gpath,"cv_pf_binary.pdf",sep=""),width=12,height=8)
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