# Set data and graphics path
dpath = "../data/"
gpath = "../graphs/"

fmri_mcmc_plot <- function(n.sim, mod, n.chains, x=1, beta=1, sigma2m=1, phi=1, sigma2s=1)
{
  # Load simulated data
  load(paste(dpath,"dlm_ar_sim-20-",mod,".rdata",sep=""))
  
  # Load MCMC data
  out.all <- list()
  for(i in 1:n.chains)
  {
    load(paste(dpath,"fmri_mcmc_test-",paste(n.sim,mod,i,beta,sigma2m,phi,sigma2s,sep="-"),".rdata",sep=""))
    out.all[[i]] = out
  }
  n.sims = out.all[[1]]$mcmc.details$n.sims
  n.burn = out.all[[1]]$mcmc.details$n.burn
  n.thin = out.all[[1]]$mcmc.details$n.thin
  n.keep = (n.sims - n.burn) %/% n.thin

  # Traceplots for beta
  if(beta)
  {
    iter = (1:n.keep)*n.thin
    d = length(mysims[[n.sim]]$true.params$beta)
    pdf(file=paste(gpath,"fmri_mcmc_test-",paste(n.sim,mod,sep="-"),"-traceplots-beta.pdf",sep=""))
    par(mfrow=c(d,1), mar=c(5,6,4,2)+.1)
    mins = apply(matrix(sapply(out.all, function(x) apply(x$beta, 2, min)),nr=d), 1, min)
    maxs = apply(matrix(sapply(out.all, function(x) apply(x$beta, 2, max)),nr=d), 1, max)

    for(i in 1:d)
    {
      ylab = bquote(expression(beta[.(i-1)]))
      plot(iter,out.all[[1]]$beta[,i],type="l",ylim=c(mins[i],maxs[i]),xlab="Iteration",ylab=eval(ylab),cex.lab=1.5)
      if(n.chains > 1)
      {
        for(j in 2:n.chains) lines(iter,out.all[[j]]$beta[,i],col=2*(j-1))
      }
      abline(h=mysims[[n.sim]]$true.params$beta[i])
    }
    dev.off()
}