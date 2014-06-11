source("pf_mc_functions.r")

# Set data and graphics path
dpath = "../data/"
gpath = "../graphs/"

# Function to plot rm pf quantiles
fmri_pl_quantiles <- function(N, mod.sim, n, nsims, nrtot, nruns, mod.est, np, alpha = 0.05, mcmc.comp, burn = 0)
{
  # Load simulated data
  load(paste(dpath,"dlm_ar_sim-",N,"-",mod.sim,"-2.rdata",sep=""))
  mysims = get(paste(mod.sim,"_dat",sep=""))[[1]][[n]]
  nt = dim(mysims[[1]]$y)[2]
  d = length(mysims[[1]]$true.params$beta)
  p.sim = length(mysims[[1]]$true.params$G[,1])
  truth.theta = c(mysims[[1]]$true.params$beta,c(mysims[[1]]$true.params$G[,1]),mysims[[1]]$true.params$W[1,1],mysims[[1]]$true.params$V[1,1])
  mod.lab = eval(bquote(expression(Sim:~M[.(paste(strsplit(mod.sim,"")[[1]][-1],sep="",collapse=""))]~Fit:~M[.(paste(strsplit(mod.est,"")[[1]][-1],sep="",collapse=""))])))

  # Load first pf run, record dimension of estimated state, and set true param labels for comparison to estimates
  load(paste(dpath,"fmri_pl-",paste(mod.sim,n,nsims[1],nrtot,mod.est,np[1],alpha,sep="-"),".rdata",sep=""))
  p = dim(pf.out$state.quant[[1]])[2]-1
  if(p.sim > p) truth = truth.theta[c(1:d,(d+1):(d+p),d+p.sim+1,d+p.sim+2)]
  if(p.sim < p) truth = c(truth.theta[c(1:d,(d+1):(d+p.sim))],rep(NA,p-p.sim),truth.theta[(d+p.sim+1):(d+p.sim+2)])
  if(p.sim == p) truth = truth.theta

  # Plot quantiles for different pf runs for each sim
  for(i in nsims)
  {
    # Find plot max and mins
    gmin = rep(Inf,d+p+3)
    gmax = rep(-Inf,d+p+3)
    for(j in 1:length(np))
    {
      load(paste(dpath,"fmri_pl-",paste(mod.sim,n,i,nrtot,mod.est,np[j],alpha,sep="-"),".rdata",sep=""))
      for(l in nruns)
      {
        if(burn > 0) tmp = pf.out$state.quant[[l]][-(1:burn),1,] else tmp = pf.out$state.quant[[l]][,1,]
        gmin[1] = min(gmin[1], tmp[,1], mysims[[i]]$x[1,],na.rm=T)
        gmax[1] = max(gmax[1], tmp[,2], mysims[[i]]$x[1,],na.rm=T)
        for(k in 2:(d+p+3))
        {
          if(burn > 0) tmp = pf.out$theta.quant[[l]][-(1:burn),k-1,] else tmp = pf.out$theta.quant[[l]][,k-1,]
          gmin[k] = min(gmin[k], tmp[,1], truth[k-1],na.rm=T)
          gmax[k] = max(gmax[k], tmp[,2], truth[k-1],na.rm=T)
        }
      }
    }
#     if(!missing(mcmc.comp))
#     {
#       mcmc.out <- list()
#       for(k in 1:mcmc.comp$n.chains)
#       {
#         load(paste(dpath,"fmri_dlm_mcmc_test-",paste(N,n,i,mod.sim,mod.est,"mle",k,mcmc.comp$nsims,mcmc.comp$nburn,mcmc.comp$nthin,"FALSE",sep="-"),".rdata",sep=""))
#         mcmc.out[[k]] = out.est
#       }
#       beta.quant = apply(sapply(mcmc.out, function(x) x$out$beta), 2, function(x) quantile(x, probs=c(.025,.5,.975)))
#       phi.quant = apply(mcmc.out2$phi[[1]], 2, function(x) quantile(x, probs=c(.025,.5,.975)))
#       sigma2s.quant = apply(mcmc.out2$sigma2s, 2, function(x) quantile(x, probs=c(.025,.5,.975)))
#       sigma2m.quant = quantile(mcmc.out2$sigma2m, probs=c(.025,.5,.975))
#       mcmc.quant2 = cbind(beta.quant,phi.quant,sigma2s.quant,sigma2m.quant)
#       x.quant2 = apply(mcmc.out2$x, 3, function(x) quantile(x, probs=c(.025,.5,.975)))
#     }
    
    # Construct plots
    # State
    pdf(file = paste(gpath,"fmri_pl-",paste(mod.sim,n,i,length(nruns),mod.est,100*alpha,sep="-"),"-states.pdf",sep=""), width = 5, height=5)
    par(mar = c(5,6,4,2) + 0.1)
    plot(0:nt, mysims[[i]]$x[1,], ylim = c(gmin[1], gmax[1]), xlab = "Time", ylab = expression(x[t]), type = "l", col = "gray", main = mod.lab, lwd = 3, cex.main = 1.5, cex.lab = 2, cex.axis=1.2)
    for(j in 1:length(np))
    {
      load(paste(dpath,"fmri_pl-",paste(mod.sim,n,i,nrtot,mod.est,np[j],alpha,sep="-"),".rdata",sep=""))
      for(k in nruns)
      {
        lines(0:nt, pf.out$state.quant[[k]][,1,1],col=j+1)
        lines(0:nt, pf.out$state.quant[[k]][,1,2],col=j+1)
      }
    }
    legend("topright",c(paste(np," particles",sep=""),"Truth"),lty=c(rep(1,length(np)),1),col=c(2:(length(np)+1),"gray"))
    dev.off()
    
    # Beta
    pdf(file = paste(gpath,"fmri_pl-",paste(mod.sim,n,i,length(nruns),mod.est,100*alpha,sep="-"),"-beta.pdf",sep=""), width = 5*d, height=5)
    par(mfrow=c(1,d),mar = c(5,6,4,2) + 0.1)
    for(k in 2:(d+1))
    {
      expr = bquote(expression(beta[.(k-2)]))
      plot(0:nt, rep(truth.theta[k-1],nt+1), ylim = c(gmin[k], gmax[k]), xlab = "", ylab = eval(expr), type = "l", col = "gray", lwd = 3, cex.lab = 2, cex.axis=1.2)
      mtext(truth.theta[k-1], at=truth.theta[k-1], side=4)
      for(j in 1:length(np))
      {
        load(paste(dpath,"fmri_pl-",paste(mod.sim,n,i,nrtot,mod.est,np[j],alpha,sep="-"),".rdata",sep=""))
        for(l in nruns)
        {
          lines(0:nt, pf.out$theta.quant[[l]][,k-1,1], col=j+1)
          lines(0:nt, pf.out$theta.quant[[l]][,k-1,2], col=j+1)
        }
      }
    }
    dev.off()
    
    # Phi
    pdf(file = paste(gpath,"fmri_pl-",paste(mod.sim,n,i,length(nruns),mod.est,100*alpha,sep="-"),"-phi.pdf",sep=""), width = 5*p, height=5)
    par(mfrow=c(1,p),mar = c(5,6,4,2) + 0.1)
    for(k in (d+2):(d+1+p))
    {
      expr = bquote(expression(phi[.(k-d-1)]))
      plot(0:nt, rep(truth.theta[k-1],nt+1), ylim = c(gmin[k], gmax[k]), xlab = "", ylab = eval(expr), type = "l", col = "gray", lwd = 3, cex.lab = 2, cex.axis=1.2)
      mtext(truth.theta[k-1], at=truth.theta[k-1], side=4)
      for(j in 1:length(np))
      {
        load(paste(dpath,"fmri_pl-",paste(mod.sim,n,i,nrtot,mod.est,np[j],alpha,sep="-"),".rdata",sep=""))
        for(l in nruns)
        {
          lines(0:nt, pf.out$theta.quant[[l]][,k-1,1], col=j+1)
          lines(0:nt, pf.out$theta.quant[[l]][,k-1,2], col=j+1)
        }
      }
    }
    dev.off()
    
    # Variances
    pdf(file = paste(gpath,"fmri_pl-",paste(mod.sim,n,i,length(nruns),mod.est,100*alpha,sep="-"),"-var.pdf",sep=""), width = 10, height=5)
    par(mfrow=c(1,2),mar = c(5,6,4,2) + 0.1)
    for(k in (d+2+p):(d+3+p))
    {
      expr = expression(sigma[s]^2,sigma[m]^2)
      plot(0:nt, rep(truth[k-1],nt+1), ylim = c(gmin[k], gmax[k]), xlab = "", ylab = expr[k-d-p-1], type = "l", col = "gray", lwd = 3, cex.lab = 2, cex.axis=1.2)
      mtext(truth[k-1], at=truth[k-1], side=4)
      for(j in 1:length(np))
      {
        load(paste(dpath,"fmri_pl-",paste(mod.sim,n,i,nrtot,mod.est,np[j],alpha,sep="-"),".rdata",sep=""))
        for(l in nruns)
        {
          lines(0:nt, pf.out$theta.quant[[l]][,k-1,1], col=j+1)
          lines(0:nt, pf.out$theta.quant[[l]][,k-1,2], col=j+1)
        }
      }
    }
    dev.off()
  }
}

require(plyr)
mydata = expand.grid(N = 20, mod.sim = c("M101","M011","M021"),n=19,nrtot=5,mod.est=c("M101","M011","M021"),stringsAsFactors = FALSE)
m_ply(mydata, function(N, mod.sim, n, nrtot, mod.est) fmri_pl_quantiles(N, mod.sim, n, 3, nrtot, 1:5, mod.est, c(100,500,1000), burn=0))

fmri_pl_loglik_wsim = function(N, mod.sim, n, nsims, nrtot, nruns, mod.est, np, alpha = 0.05)
{
  # Plot marginal likelihoods for each sim for each # particles
  for(n.sim in nsims)
  {
    # Load simulated data
    load(paste(dpath,"dlm_ar_sim-",N,"-",mod.sim,"-2.rdata",sep=""))
    mysims = get(paste(mod.sim,"_dat",sep=""))[[1]][[n]]
    d = length(mysims[[n.sim]]$true.params$beta)
    p = length(mysims[[n.sim]]$true.params$G[,1])
    truth.theta = c(mysims[[n.sim]]$true.params$beta,mysims[[n.sim]]$true.params$G[,1],mysims[[n.sim]]$true.params$W[1,1],mysims[[n.sim]]$true.params$V[1,1])

    # Parameter and model labels
    mlabels = rep(NA,length(mod.est))
    for(i in 1:length(mod.est)) mlabels[i] = eval(bquote(expression(M[.(paste(strsplit(mod.est[i],"")[[1]][-1],sep="",collapse=""))])))
    blabels = rep(NA, d+p+2)
    for(i in 1:d) blabels[i] = eval(bquote(expression(paste(beta[.(i-1)]," = ",.(truth.theta[i]),sep = ""))))
    for(i in 1:p) blabels[i+d] = eval(bquote(expression(paste(phi[.(i)]," = ",.(truth.theta[i+d]),sep=""))))
    blabels[d+p+1] = eval(bquote(expression(paste(sigma[s]^2," = ",.(truth.theta[d+p+1]),sep=""))))
    blabels[d+p+2] = eval(bquote(expression(paste(sigma[m]^2," = ",.(truth.theta[d+p+2]),sep=""))))
    
    # Load approximate log marginal likelihoods for particle filter runs
    lmargliks = list()
    for(j in 1:length(np))
    {
      lmargliks[[j]] = matrix(nr = nruns, nc = length(mod.est))
      for(k in 1:length(mod.est))
      {
        load(paste(dpath,"fmri_pl-",paste(mod.sim,n,n.sim,nrtot,mod.est[k],np[j],alpha,sep="-"),".rdata",sep=""))
        lmargliks[[j]][,k] = pf.out$lmarglik
      }
    }
    
    # Plot kernel density estimates of of log marginal likelihoods for pfs compared with truth
    dmax = max(sapply(lmargliks, function(x) apply(x, 2, function(a) max(density(a)$y,na.rm=TRUE))))
    bmax = max(sapply(lmargliks, function(x) apply(x, 2, function(a) max(density(a)$x,na.rm=TRUE))))   
    bmin = min(sapply(lmargliks, function(x) apply(x, 2, function(a) min(density(a)$x,na.rm=TRUE))))
    cols = rainbow(length(np))
    pdf(file=paste(gpath,"fmri_pl_loglik-",paste(mod.sim,n,n.sim,nruns,sep="-"),".pdf",sep=""),width=5*length(mod.est),height=5)
    par(mfrow=c(1,length(mod.est)),mar=c(5,6,4,2)+0.1)
    xlab = c("Log marginal likelihood",rep("",length(mod.est)-1))
    ylab = c("Density",rep("",length(mod.est)-1))
    for(i in 1:length(mod.est))
    {
      plot(density(lmargliks[[1]][,i]),lwd=2,col=cols[1],main=mlabels[i],xlab=xlab[i],ylab=ylab[i],xlim=c(bmin,bmax),ylim=c(0,dmax),cex.axis=1.5,cex.lab=1.75,cex.main=2)
      if(length(np) > 1) for(k in 2:length(np)) lines(density(lmargliks[[k]][,i]),lwd=2,col=cols[k])
      if(i == 1) legend("topleft",legend=paste(np," particles",sep=""),lty=rep(1,length(np)),lwd=rep(2,length(np)),col=cols,cex=1.5)
      if(i == which(mod.est == mod.sim)) for(j in 1:(d+p+2)) mtext(blabels[j],side=3,at=bmin + ((j-1)/(d+p+2))*(bmax-bmin),cex=0.85)
    }
    dev.off()
  }
}

require(plyr)
mydata = expand.grid(N = 20, mod.sim = "M101", n = 19, nsims = 3, nrtot = 5, nruns = 5)
m_ply(mydata, function(N, mod.sim, n, nsims, nrtot, nruns) fmri_pl_loglik_wsim(N, mod.sim, n, nsims, nrtot, nruns, c("M101","M011","M021"), c(100,500,1000)))

fmri_pl_comp_wsim = function(N, mod.sim, n, nsims, nrtot, nruns, mod.est, np, alpha = 0.05)
{
  require(compositions)
  
  # Plot quantiles for each sim for each # particles
  for(n.sim in nsims)
  {
    # Load simulated data
    load(paste(dpath,"dlm_ar_sim-",N,"-",mod.sim,"-2.rdata",sep=""))
    mysims = get(paste(mod.sim,"_dat",sep=""))[[1]][[n]]
    truth.theta = c(mysims[[n.sim]]$true.params$beta,mysims[[n.sim]]$true.params$G[,1],mysims[[n.sim]]$true.params$W[1,1],mysims[[n.sim]]$true.params$V[1,1])
    
    # Load approximate log marginal likelihoods for particle filter runs
    lmargliks = list()
    for(j in 1:length(np))
    {
      lmargliks[[j]] = matrix(nr = nruns, nc = length(mod.est))
      for(k in 1:length(mod.est))
      {
        load(paste(dpath,"fmri_pl-",paste(mod.sim,n,n.sim,nrtot,mod.est[k],np[j],alpha,sep="-"),".rdata",sep=""))
        lmargliks[[j]][,k] = pf.out$lmarglik
      }
    }
    
    # Ternary diagrams of posterior model probabilities
    # Compute pl approximate posterior model probabilities
    pmargliks = list()
    for(i in 1:length(np))
    {
      pmargliks[[i]] = apply(lmargliks[[i]], 1, function(a) postModProbs(a, rep(1/length(mod.est), length(mod.est))))
      pmargliks[[i]] = t(pmargliks[[i]])
    }
    
    # Ternary diagrams of posterior model probabilities
    require(compositions)
    mlabels = rep(NA,length(mod.est))
    for(i in 1:length(mod.est)) mlabels[i] = eval(bquote(expression(M[.(paste(strsplit(mod.est[i],"")[[1]][-1],sep="",collapse=""))])))
    pdf(file=paste(gpath,"fmri_pl_ternary-",paste(mod.sim,n,n.sim,nruns,sep="-"),".pdf",sep=""),width=10,height=10)
    par(mfrow=c(2,2))
    for(i in 1:length(np))
    {
      plot(acomp(pmargliks[[i]]),labels=mlabels, lwd=2)
      mtext(paste(np[i]," particles",sep=""),side=3,cex=2)
    }
    dev.off()
  }
}

require(plyr)
mydata = expand.grid(N = 20, mod.sim = c("M101","M011","M021"), n = 19, nsims = 3, nrtot = 5, nruns = 5)
m_ply(mydata, function(N, mod.sim, n, nsims, nrtot, nruns) fmri_pl_comp_wsim(N, mod.sim, n, nsims, nrtot, nruns, c("M101","M011","M021"), c(100,500,1000)))

