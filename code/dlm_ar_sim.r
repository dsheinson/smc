source("dlm_sim.r")
source("dlm_ar_functions.r")

# Set data and graphics path
dpath = "../data/"
gpath = "../graphs/"

# Simulate data from dynamic regression model
set.seed(61)
dlm_ar_sim <- function(N, mod, beta, sigma2m = 0, phi = NULL, sigma2s = 0, rho = NULL, sigma2b = 0, label="")
{
  conv = scan(paste(dpath,"basis_sim-245.txt",sep=""))
  nt = length(conv)
  U = array(rbind(1,conv), c(1,2,nt))
  beta = rep(beta,1); stopifnot(length(beta) == 2)

  if(mod == "dr")
  {
    p = length(rho)
    F = array(0, c(1,p,nt)); F[1,1,] = conv
    V = matrix(sigma2m, nr = 1, nc = 1)
    G = makeG(rho)
    W = sigma2b*(c(1,rep(0,p-1))%*%t(c(1,rep(0,p-1))))
    C0 = makeC0(rho)
    x0 = rep(t(chol(C0))%*%rnorm(p,0,1),1)
  } else if(mod == "ar") {
    p = length(phi)
    F = array(0, c(1,p,nt)); F[1,1,] = 1
    V = matrix(sigma2m, nr = 1, nc = 1)
    G = makeG(phi)
    W = sigma2s*(c(1,rep(0,p-1))%*%t(c(1,rep(0,p-1))))
    C0 = makeC0(phi)
    x0 = rep(t(chol(C0))%*%rnorm(p,0,1),1)
  } else if(mod == "both") {
    p1 = length(rho); p2 = length(phi)
    F = array(0, c(1,p1+p2,nt)); F[1,1,] = conv; F[1,1+p1,] = 1
    V = matrix(sigma2m, nr = 1, nc = 1)
    G = bdiag(list(makeG(rho),makeG(phi)))
    eb = c(sigma2b,rep(0,p1-1))
    es = c(sigma2s,rep(0,p2-1))
    W = bdiag(list(eb%*%t(eb),es%*%t(es)))
    C0 = bdiag(list(makeC0(rho),makeC0(phi)))
    x0 = rep(t(chol(C0))%*%rnorm(p1+p2,0,1),1)
  } else { stop("mod must be 'dr','ar', or 'both'")}

  mysims = list()
  for(j in 1:N) mysims[[j]] = dlm.sim(nt, F, G, V, W, x0, beta, U)

  # Save data
  save(mysims, file=paste(dpath,"dlm_ar_sim-",N,"-",label,".rdata",sep=""))
}

dlm_ar_sim(20, "dr", label="M101", beta = c(900,1), sigma2m = 1, rho = .9, sigma2b = 1)
dlm_ar_sim(20, "ar", label="M010", beta = c(900,1), phi = .5, sigma2s = 1)
dlm_ar_sim(20, "ar", label="M011", beta = c(900,1), sigma2m = 1, phi = .5, sigma2s = 1)
dlm_ar_sim(20, "ar", label="M020", beta = c(900,1), phi = c(.5,.25), sigma2s = 1)
dlm_ar_sim(20, "both", label="M111", beta = c(900,1), sigma2m = 1, phi = .5, sigma2s = 1, rho = .9, sigma2b = 1)

# Plot data
dlm_ar_plot <- function(N, label)
{
  load(paste(dpath,"dlm_ar_sim-",N,"-",label,".rdata",sep=""))
  if(label == "M101") title = expression(paste("Simulated time series from ",M[101],sep=""))
  if(label == "M010") title = expression(paste("Simulated time series from ",M[010],sep=""))
  if(label == "M011") title = expression(paste("Simulated time series from ",M[011],sep=""))
  if(label == "M020") title = expression(paste("Simulated time series from ",M[020],sep=""))
  if(label == "M111") title = expression(paste("Simulated time series from ",M[111],sep=""))
  
  for(i in 1:N)
  {
    nt = dim(mysims[[i]]$y)[2]
    pdf(file=paste(gpath,"dlm_ar_sim-",i,"-",label,".pdf",sep=""))
    par(mfrow=c(2,1),mar=c(5,5,4,2)+0.1)
    plot(0:nt,c(NA,mysims[[i]]$y[1,]),type="l",xlab="",ylab=expression(y[t]),cex.lab=1.5)
    title(title)
    if(label == "M111"){
      ymin = min(mysims[[i]]$x); ymax = max(mysims[[i]]$x)
      plot(0:nt,mysims[[i]]$x[1,],type="l",ylim = c(ymin,ymax),xlab=expression(t),ylab=expression(x[t]),cex.lab=1.5)     
      lines(0:nt,mysims[[i]]$x[2,],col=2)
    } else {
      plot(0:nt,mysims[[i]]$x[1,],type="l",xlab=expression(t),ylab=expression(x[t]),cex.lab=1.5)
    }
    dev.off()
  }
}

require(plyr)
mydata = data.frame(N=rep(20,5),label=c("M101","M010","M011","M020","M111"),stringsAsFactors=FALSE)
m_ply(mydata, dlm_ar_plot)
