source("pf_mc_functions.r")
source("rm_fmri_functions.r")
require(plyr)

# Set data and graphics path
dpath = "../data/"
gpath = "../graphs/"

# Function to plot rm pf quantiles
fmri_rm_quantiles <- function(N, mod.sim, dimx, n, nsims, nruns, mod.est, np, sd.fac = 1, alpha = 0.05, burn = 10)
{
  # Load simulated data
  load(paste(dpath,"dlm_ar_sim-",N,"-",mod.sim,"-",dimx,".rdata",sep=""))
  mysims = get(paste(mod.sim,"_dat",sep=""))[[1]][[n]]
  nt = dim(mysims[[1]]$y)[2]
  expr = expression(beta[0],beta[1],phi,sigma[s]^2,sigma[m]^2)
  
  # Plot quantiles for different rm runs for each sim, compared with MCMC and MLE
  gmin = rep(Inf,6)
  gmax = rep(-Inf,6)
  for(j in 1:nsims)
  {
    # Load MCMC and MLEs
    out.all = list()
    for(k in 1:3)
    {
      load(paste(dpath,"fmri_dlm_mcmc_test-",paste(N,n,j,mod.sim,dimx,k,11000,1000,10,FALSE,sep="-"),".rdata",sep=""))
      out.all[[k]] = out.est
    }
    x.mcmc <- quantile(sapply(out.all, function(a) a$out$x[,1,nt+1]), c(alpha/2,1-alpha/2))
    beta0.mcmc <- quantile(sapply(out.all, function(x) x$out$beta[,1]), c(alpha/2,1-alpha/2))
    beta1.mcmc <- quantile(sapply(out.all, function(x) x$out$beta[,2]), c(alpha/2,1-alpha/2))
    phi.mcmc <- quantile(sapply(out.all, function(x) x$out$phi[[1]][,1]), c(alpha/2,1-alpha/2))
    sigma2m.mcmc <- quantile(sapply(out.all, function(x) x$out$sigma2m), c(alpha/2,1-alpha/2))
    sigma2s.mcmc <- quantile(sapply(out.all, function(x) x$out$sigma2s[,1]), c(alpha/2,1-alpha/2))
    mcmc.quant <- rbind(x.mcmc,beta0.mcmc,beta1.mcmc,phi.mcmc,sigma2s.mcmc,sigma2m.mcmc)
    truth.theta = c(mysims[[j]]$true.params$beta,mysims[[j]]$true.params$G[1,1],mysims[[j]]$true.params$W[1,1],mysims[[j]]$true.params$V[1,1])
    
    # Find plot max and mins
    for(i in 1:nruns)
    {
      load(paste(dpath,"fmri_rm-",paste(N,mod.sim,dimx,n,j,i,mod.est,np,sd.fac,100*alpha,sep="-"),".rdata",sep=""))
      gmin[1] = min(gmin[1], pf.out$state.quant[-(1:burn),1,1], mysims[[j]]$x[1,],mcmc.quant[1,1],out.all[[1]]$x.mle[1,])
      gmax[1] = max(gmax[1], pf.out$state.quant[-(1:burn),1,2], mysims[[j]]$x[1,],mcmc.quant[1,2],out.all[[1]]$x.mle[1,])
      for(k in 2:6)
      {
        gmin[k] = min(gmin[k], pf.out$theta.quant[-(1:burn),k-1,1], mcmc.quant[k,1], truth.theta[k-1], out.all[[1]]$theta.mle[k-1])
        gmax[k] = max(gmax[k], pf.out$theta.quant[-(1:burn),k-1,2], mcmc.quant[k,2], truth.theta[k-1], out.all[[1]]$theta.mle[k-1])
      }
      rm(pf.out)
    }
    
    # Construct plots
    pdf(file =paste(gpath,"fmri_rm-",paste(N,mod.sim,dimx,n,j,nruns,mod.est,np,100*sd.fac,100*alpha,sep="-"),".pdf",sep=""), width = 15, height=10)
    par(mfrow = c(2,3), mar = c(5,6,4,2) + 0.1)
    plot(0:nt, out.all[[1]]$x.mle[1,], ylim = c(gmin[1], gmax[1]), xlab = "Time", ylab = expression(x[t]), type = "l", col = "gray", lwd = 3, cex.lab = 1.5)
    for(i in 1:nruns)
    {
      load(paste(dpath,"fmri_rm-",paste(N,mod.sim,dimx,n,j,i,mod.est,np,sd.fac,100*alpha,sep="-"),".rdata",sep=""))
      lines(0:nt, pf.out$state.quant[,1,1], col = 2)
      lines(0:nt, pf.out$state.quant[,1,2], col = 2)
      rm(pf.out)
    }
    lines(0:nt, mysims[[j]]$x[1,])
    points(c(nt,nt),mcmc.quant[1,],col=4,pch="x",cex=2)
    legend("topright",c("95% CI","MLE","Truth","MCMC"),lty=c(1,1,1,NA),pch = c(NA,NA,NA,"x"),col=c(2,"gray",1,4))
    for(k in 2:6)
    {
      plot(0:nt, rep(out.all[[1]]$theta.mle[k-1],nt+1), ylim = c(gmin[k], gmax[k]), xlab = "", ylab = expr[k-1], type = "l", col = "gray", lwd = 3, cex.lab = 1.5)
      for(i in 1:nruns)
      {
        load(paste(dpath,"fmri_rm-",paste(N,mod.sim,dimx,n,j,i,mod.est,np,sd.fac,100*alpha,sep="-"),".rdata",sep=""))
        lines(0:nt, pf.out$theta.quant[,k-1,1], col = 2)
        lines(0:nt, pf.out$theta.quant[,k-1,2], col = 2)
      }
      abline(h=truth.theta[k-1])
      points(c(nt,nt),mcmc.quant[k,],col=4,pch="x",cex=2)
    }
    dev.off()
  }
}

mydata = expand.grid(N = 20, mod.sim = c("M101","M011"),dimx=2,n=6,nsims=8,nruns=3,mod.est=c("M101","M011"),np=100,sd.fac=1,alpha=0.05,burn=10,stringsAsFactors = FALSE)
m_ply(mydata, fmri_rm_quantiles)

# Function to plot rm pf posterior model probabilities
fmri_rm_lik <- function(N, mod.sim, dimx, n, nsims, nruns, np, sd.fac, byn, ylim, prior = FALSE, alpha = 0.05)
{
  # Load simulated data
  load(paste(dpath,"dlm_ar_sim-",N,"-",mod.sim,"-",dimx,".rdata",sep=""))
  mysims = get(paste(mod.sim,"_dat",sep=""))[[1]]
 
  # Load approximate log marginal likelihoods for particle filter runs
  rm.lmarglik = array(NA, c(nruns, nsims, 2, length(n), length(sd.fac)))
  for(l in 1:length(n))
  {
    for(m in 1:length(sd.fac))
    {
      for(i in 1:nruns)
      {
        for(j in 1:nsims)
        {
          for(k in 1:2)
          {
            mod = c("M101","M011")[k]
            if(!prior | m == 1)
            {
              load(paste(dpath,"fmri_rm-",paste(N,mod.sim,dimx,n[l],j,i,mod,np,sd.fac[m],100*alpha,sep="-"),".rdata",sep=""))
              rm.lmarglik[i,j,k,l,m] = pf.out$lmarglik
            } else {
              require(smcUtils)
              load(paste(dpath,"fmri_rm-",paste(N,mod.sim,dimx,n[l],j,i,mod,np,sd.fac[1],100*alpha,sep="-"),".rdata",sep=""))
              prior.old = pf.out$prior
              prior.new = rprior.convert(pf.out$wt.mom$center, diag(pf.out$wt.mom$cov), sd.fac[m])
              dlprior.old = function(x, theta) dlprior(x, theta, prior.old)
              dlprior.new = function(x, theta) dlprior(x, theta, prior.new)
              rm.lmarglik[i,j,k,l,m] = pf.lmarglik.prior(pf.out$out, dlprior.old, dlprior.new)
              print(paste(i,j,k,l,m,sep="-"))
            }
          }
        }
      }
    }
  }
   
  # Input y-limits?
  if(missing(ylim)) ylim = c(min(rm.lmarglik),max(rm.lmarglik))
  
  # Plot log marginal likelihoods estimated by M101 vs M011
  if(!prior) file=paste(gpath,"fmri_rm_lik-",paste(N, mod.sim, dimx, nsims, nruns, np, byn, sep="-"),".pdf",sep="") else file=paste(gpath,"fmri_rm_lik-",paste(N, mod.sim, dimx, nsims, nruns, np, byn, "prior", sep="-"),".pdf",sep="")
  pdf(file=file)
  par(mfrow=c(2,2),mar=c(5,6,4,2)+0.1)
  if(byn)
  {
    for(i in 1:length(n))
    {
      beta = mysims[[n[i]]][[1]]$true.params$beta
      phi = mysims[[n[i]]][[1]]$true.params$G[1,1]
      sigma2s = mysims[[n[i]]][[1]]$true.params$W[1,1]
      sigma2m = mysims[[n[i]]][[1]]$true.params$V[1,1]
      if(i == 1)
      {
        xlab = "Simulation"
        ylab = "Log marginal likelihood"
        main = substitute(paste(aa," simulation, ",ab," particles",sep=""),list(aa=mod.sim,ab=np))
        mtext = substitute(paste(beta[0]," = ",aa,", ",beta[1]," = ",ab,", ",phi," = ",ac,", ",sigma[s]^2," = ",ad,", ",sigma[m]^2," = ",ae,sep=""),list(aa=beta[1],ab=beta[2],ac=phi,ad=sigma2s,ae=sigma2m))
      } else {
        xlab = ylab = main = ""
        mtext = substitute(paste(beta[0]," = ",aa,", ",beta[1]," = ",ab,", ",phi," = ",ac,", ",sigma[s]^2," = ",ad,", ",sigma[m]^2," = ",ae,sep=""),list(aa=beta[1],ab=beta[2],ac=phi,ad=sigma2s,ae=sigma2m))
      }
      plot(rep(1,nruns), rm.lmarglik[,1,1,i,1], xlim = c(1,nsims), ylim = ylim, col = 2, xlab = xlab, ylab = ylab, main = main, cex.lab = 1.5)
      mtext(mtext,side=3,cex=0.85)
      if(nsims > 1) for(j in 2:nsims) points(rep(j,nruns), rm.lmarglik[,j,1,i,1], col=2)
      for(j in 1:nsims) points(rep(j,nruns), rm.lmarglik[,j,2,i,1], col = 4)
      if(i == 1) legend("bottomright",c("M101","M011"),pch=c(1,1),col=c(2,4))
    }
  } else {
    beta = mysims[[n[1]]][[1]]$true.params$beta
    phi = mysims[[n[1]]][[1]]$true.params$G[1,1]
    sigma2s = mysims[[n[1]]][[1]]$true.params$W[1,1]
    sigma2m = mysims[[n[1]]][[1]]$true.params$V[1,1]
    for(i in 1:length(sd.fac))
    {
      ind.M101 = (mod.sim != "M101") + i*(mod.sim == "M101")
      ind.M011 = (mod.sim != "M011") + i*(mod.sim == "M011")
      if(i == 1)
      {
        xlab = "Simulation"
        ylab = "Log marginal likelihood"
        main = substitute(paste(aa," simulation, ",ab," particles",sep=""),list(aa=mod.sim,ab=np))
        mtext = substitute(paste(beta[0]," = ",aa,", ",beta[1]," = ",ab,", ",phi," = ",ac,", ",sigma[s]^2," = ",ad,", ",sigma[m]^2," = ",ae,sep=""),list(aa=beta[1],ab=beta[2],ac=phi,ad=sigma2s,ae=sigma2m))
      } else {
        xlab = ylab = ""
        main = paste(mod.sim," SD Factor = ",sd.fac[i],sep="")
        mtext = substitute(paste(beta[0]," = ",aa,", ",beta[1]," = ",ab,", ",phi," = ",ac,", ",sigma[s]^2," = ",ad,", ",sigma[m]^2," = ",ae,sep=""),list(aa=beta[1],ab=beta[2],ac=phi,ad=sigma2s,ae=sigma2m))
      }
      plot(rep(1,nruns), rm.lmarglik[,1,1,1,ind.M101], xlim = c(1,nsims), ylim = ylim, col = 2, xlab = xlab, ylab = ylab, main = main, cex.lab = 1.5)
      mtext(mtext,side=3,cex=0.85)
      if(nsims > 1) for(j in 2:nsims) points(rep(j,nruns), rm.lmarglik[,j,1,1,ind.M101], col=2)
      for(j in 1:nsims) points(rep(j,nruns), rm.lmarglik[,j,2,1,ind.M011], col = 4)
      if(i == 1) legend("bottomright",c("M101","M011"),pch=c(1,1),col=c(2,4))
    }
  }
  dev.off()
}

mydata = expand.grid(N = 20, mod.sim = c("M101","M011"),dimx=2,nsims=8,nruns=3,np=100,stringsAsFactors = FALSE)
m_ply(mydata, function(N,mod.sim,dimx,nsims,nruns,np) fmri_rm_lik(N, mod.sim, dimx, n=6, nsims, nruns, np, sd.fac=c(1,3,5,7), byn=FALSE, alpha=0.05))

mydata = expand.grid(N = 20, mod.sim = c("M101","M011"),dimx=2,nsims=8,nruns=3,np=100,stringsAsFactors = FALSE)
m_ply(mydata, function(N,mod.sim,dimx,nsims,nruns,np) fmri_rm_lik(N, mod.sim, dimx, n=6, nsims, nruns, np, sd.fac=c(1,3,5,7), byn=FALSE, prior = TRUE, alpha=0.05))

mydata = expand.grid(N = 20, mod.sim = c("M101","M011"), dimx=2, nsims=8, nruns=3, np=100, stringsAsFactors = FALSE)
m_ply(mydata, function(N,mod.sim,dimx,nsims,nruns,np) fmri_rm_lik(N, mod.sim, dimx, n=c(1,6,11), nsims, nruns, np, sd.fac=1, byn=TRUE, alpha=0.05))
