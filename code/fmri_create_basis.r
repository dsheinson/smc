source("fmri_sim.r")

# Set data and graphics path
dpath = "../data/"
gpath = "../graphs/"

# Function to simulate design matrices and choose most stable
optim.design <- function(N, n, E, intercept = TRUE)
{
  hrf = boynton.hrf(n)$hrf
  designs = list(); length(designs) = N
  for(j in 1:N)
  {
    boxcar = rapidEvent.boxcar(n, E)
    basis1 = fmri.convolve(hrf,boxcar$boxcar[,1])
    if(E > 1)
    {
      basis = matrix(NA, nr = length(basis1$t), nc = E)
      basis[,1] = basis1$basis
      for(i in 2:E) basis[,i] = fmri.convolve(hrf,boxcar$boxcar[,i])$basis
    } else {basis = as.matrix(basis1)}
    if(intercept) X = cbind(1,basis) else X = basis
    eig = eigen(t(X)%*%X)
    designs[[j]] = list(boxcar=boxcar, X=X, min.eig = min(eig$values))
    cat(j,eig$values,"\n")
  }
  
  # Find optimal design
  ind = which(sapply(designs, function(x) x$min.eig) == max(sapply(designs, function(x) x$min.eig)))
  
  # Plot optimal design
  pdf(file = paste(gpath,"fmri-optim-design.pdf",sep=""))
  par(mfrow=c(3,1))
  ymax = max(designs[[ind]]$boxcar$boxcar)
  plot(designs[[ind]]$boxcar$t,designs[[ind]]$boxcar$boxcar[,1],ylim=c(0,ymax),type="l",xlab="",ylab="Boxcar")
  if(E > 1) for(i in 2:E) lines(designs[[ind]]$boxcar$t,designs[[ind]]$boxcar$boxcar[,i],col=i)
  plot(designs[[ind]]$boxcar$t,hrf,type="l",xlab="",ylab="HRF")
  plot(basis1$t,designs[[ind]]$X[,intercept+1],type="l",ylim=c(min(designs[[ind]]$X),max(designs[[ind]]$X)),xlab="Convolution",ylab="Time (s)")
  if(E > 1) for(i in 2:E) lines(basis1$t,designs[[ind]]$X[,intercept+i],col=i)
  dev.off()
}

# # Generate basis function
# set.seed(61)
# n = 245
# hrf.raw = read.table(paste(dpath,"hrf.csv",sep=""),sep=",",header=FALSE)[,1]
# boxcar = rapidEvent.boxcar(n,10)
# hrf = c(hrf.raw,rep(0,n - length(hrf.raw)))
# basis = fmri.convolve(hrf,boxcar)
# 
# # Save basis function
# write(basis,file=paste(dpath,"basis_sim-",n,".txt",sep=""))
