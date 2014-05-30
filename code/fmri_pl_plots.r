source("pf_mc_functions.r")

# Set data and graphics path
dpath = "../data/"
gpath = "../graphs/"

# Function to plot rm pf quantiles
fmri_pl_quantiles <- function(N, mod.sim, n, nsims, nrtot, nruns, mod.est, np, alpha = 0.05, burn = 0)
{
  # Load simulated data
  load(paste(dpath,"dlm_ar_sim-",N,"-",mod.sim,"-2.rdata",sep=""))
  mysims = get(paste(mod.sim,"_dat",sep=""))[[1]][[n]]
  nt = dim(mysims[[1]]$y)[2]
  d = length(mysims[[1]]$true.params$beta)

  # Plot quantiles for different pf runs for each sim
  for(j in 1:nsims)
  {
    load(paste(dpath,"fmri_pl-",paste(N,mod.sim,n,j,nrtot,mod.est,np,alpha,sep="-"),".rdata",sep=""))
    truth.theta = c(mysims[[j]]$true.params$beta,mysims[[j]]$true.params$G[1,1],mysims[[j]]$true.params$W[1,1],mysims[[j]]$true.params$V[1,1])
    
    # Find plot max and mins
    gmin = rep(Inf,d+4)
    gmax = rep(-Inf,d+4)
    for(i in 1:nruns)
    {
      if(burn > 0) tmp = pf.out$state.quant[[i]][-(1:burn),1,] else tmp = pf.out$state.quant[[i]][,1,]
      gmin[1] = min(gmin[1], tmp[,1], mysims[[j]]$x[1,])
      gmax[1] = max(gmax[1], tmp[,2], mysims[[j]]$x[1,])
      for(k in 2:(d+4))
      {
        if(burn > 0) tmp = pf.out$theta.quant[[i]][-(1:burn),k-1,] else tmp = pf.out$theta.quant[[i]][,k-1,]
        gmin[k] = min(gmin[k], tmp[,1], truth.theta[k-1])
        gmax[k] = max(gmax[k], tmp[,2], truth.theta[k-1])
      }
    }
    
    # Construct plots
    pdf(file = paste(gpath,"fmri_rm-","fmri_pl-",paste(N,mod.sim,n,j,nrtot,nruns,mod.est,np,100*alpha,sep="-"),".pdf",sep=""), width = 15, height=10)
    par(mfrow = c(2,3), mar = c(5,6,4,2) + 0.1)
    main = paste('Sim: ', mod.sim,", Fit: ",mod.est,", Particles: ",np,sep="")
    plot(0:nt, mysims[[j]]$x[1,], ylim = c(gmin[1], gmax[1]), xlab = "Time", ylab = expression(x[t]), type = "l", col = "gray", main = main, lwd = 3, cex.main = 2, cex.lab = 2, cex.axis=1.2)
    for(i in 1:nruns)
    {
      lines(0:nt, pf.out$state.quant[[i]][,1,1])
      lines(0:nt, pf.out$state.quant[[i]][,1,2])
    }
    legend("topright",c("95% CI","Truth"),lty=c(1,1),col=c("gray",1),cex=1.3)
    for(k in 2:(d+1))
    {
      expr = bquote(expression(beta[.(k-2)]))
      plot(0:nt, rep(truth.theta[k-1],nt+1), ylim = c(gmin[k], gmax[k]), xlab = "", ylab = eval(expr), type = "l", col = "gray", lwd = 3, cex.lab = 2, cex.axis=1.2)
      for(i in 1:nruns)
      {
        lines(0:nt, pf.out$theta.quant[[i]][,k-1,1])
        lines(0:nt, pf.out$theta.quant[[i]][,k-1,2])
      }
    }
    expr = expression(phi,sigma[s]^2,sigma[m]^2)
    for(k in (d+2):(d+4))
    {
      plot(0:nt, rep(truth.theta[k-1],nt+1), ylim = c(gmin[k], gmax[k]), xlab = "", ylab = expr[k-d-1], type = "l", col = "gray", lwd = 3, cex.lab = 2, cex.axis=1.2)
      for(i in 1:nruns)
      {
        lines(0:nt, pf.out$theta.quant[[i]][,k-1,1])
        lines(0:nt, pf.out$theta.quant[[i]][,k-1,2])
      }
    }
    dev.off()
  }
}

require(plyr)
mydata = expand.grid(N = 20, mod.sim = "M011",n=15,nsims=1,nrtot=20,nruns=5,mod.est=c("M101","M011","M101s"),np=c(100,500,1000,5000),alpha=0.05,burn=0,stringsAsFactors = FALSE)
m_ply(mydata, fmri_pl_quantiles)

fmri_pl_loglik_wsim = function(N, mod.sim, n, nsims, nrtot, nruns, mod.est, np, alpha = 0.05)
{
  # Plot marginal likelihoods for each sim for each # particles
  for(n.sim in 1:nsims)
  {
    # Load simulated data
    load(paste(dpath,"dlm_ar_sim-",N,"-",mod.sim,"-2.rdata",sep=""))
    mysims = get(paste(mod.sim,"_dat",sep=""))[[1]][[n]]

    # Parameter and model labels
    mlabels = rep(NA,length(mod.est))
    for(i in 1:length(mod.est)) mlabels[i] = eval(bquote(expression(M[.(paste(strsplit(mod.est[i],"")[[1]][-1],sep="",collapse=""))])))
    truth.theta = c(mysims[[n.sim]]$true.params$beta,mysims[[n.sim]]$true.params$G[1,1],mysims[[n.sim]]$true.params$W[1,1],mysims[[n.sim]]$true.params$V[1,1])
 
    # Load approximate log marginal likelihoods for particle filter runs
    lmargliks = list()
    for(j in 1:length(np))
    {
      lmargliks[[j]] = matrix(nr = nruns, nc = length(mod.est))
      for(k in 1:length(mod.est))
      {
        load(paste(dpath,"fmri_pl-",paste(N,mod.sim,n,n.sim,nrtot,mod.est[k],np[j],alpha,sep="-"),".rdata",sep=""))
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
      if(length(np > 1)) for(k in 2:length(np)) lines(density(lmargliks[[k]][,i]),lwd=2,col=cols[k])
      if(i == 1) legend("topleft",legend=paste(np," particles",sep=""),lty=rep(1,length(np)),lwd=rep(2,length(np)),col=cols,cex=1.5)
      if(i == which(mod.est == mod.sim)) mtext(substitute(paste(beta[0]," = ",aa,", ",beta[1]," = ",ab,", ",phi," = ",ac,", ",sigma[s]^2," = ",ad,", ",sigma[m]^2," = ",ae,sep=""),list(aa=truth.theta[1],ab=truth.theta[2],ac=truth.theta[3],ad=truth.theta[4],ae=truth.theta[5])),side=3,cex=0.85)
    }
    dev.off()
  }
}

require(plyr)
mydata = expand.grid(N = 20, mod.sim = "M011", n = 15, nsims = 1, nrtot = 20, nruns = 20)
m_ply(mydata, function(N, mod.sim, n, nsims, nrtot, nruns) fmri_pl_loglik_wsim(N, mod.sim, n, nsims, nrtot, nruns, c("M101","M011","M101s"), c(100, 500, 1000,5000)))

fmri_pl_comp_wsim = function(N, mod.sim, n, nsims, nrtot, nruns, mod.est, np, alpha = 0.05)
{
  require(compositions)
  
  # Plot quantiles for each sim for each # particles
  for(n.sim in 1:nsims)
  {
    # Load simulated data
    load(paste(dpath,"dlm_ar_sim-",N,"-",mod.sim,"-2.rdata",sep=""))
    mysims = get(paste(mod.sim,"_dat",sep=""))[[1]][[n]]
    truth.theta = c(mysims[[n.sim]]$true.params$beta,mysims[[n.sim]]$true.params$G[1,1],mysims[[n.sim]]$true.params$W[1,1],mysims[[n.sim]]$true.params$V[1,1])
    
    # Load approximate log marginal likelihoods for particle filter runs
    lmargliks = list()
    for(j in 1:length(np))
    {
      lmargliks[[j]] = matrix(nr = nruns, nc = length(mod.est))
      for(k in 1:length(mod.est))
      {
        load(paste(dpath,"fmri_pl-",paste(N,mod.sim,n,n.sim,nrtot,mod.est[k],np[j],alpha,sep="-"),".rdata",sep=""))
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
mydata = expand.grid(N = 20, mod.sim = "M011", n = 15, nsims = 1, nrtot = 20, nruns = 20)
m_ply(mydata, function(N, mod.sim, n, nsims, nrtot, nruns) fmri_pl_comp_wsim(N, mod.sim, n, nsims, nrtot, nruns, c("M101","M011","M101s"), c(100, 500, 1000, 5000)))

