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
mydata = expand.grid(N = 20, mod.sim = "M011",n=15,nsims=1,nrtot=20,nruns=5,mod.est=c("M101","M011","M101s"),np=c(100,500,1000),alpha=0.05,burn=0,stringsAsFactors = FALSE)
m_ply(mydata, fmri_pl_quantiles)
