source("dlm_sim.r")
source("dlm_ar_functions.r")

# Set data and graphics path
dpath = "../data/"
gpath = "../graphs/"

# Simulate data from dynamic regression model
#set.seed(61)
dlm_ar_sim <- function(N, mod, beta, sigma2m = 0, phi = NULL, sigma2b = 0, rho = NULL, sigma2s = 0)
{
  conv = scan(paste(dpath,"basis_sim-245.txt",sep=""))
  nt = length(conv)
  U = array(rbind(1,conv), c(1,2,nt))
  beta = rep(beta,1); stopifnot(length(beta) == 2)

  if(mod == "dr")
  {
    p = length(phi)
    F = array(0, c(1,p,nt)); F[1,1,] = conv
    V = matrix(sigma2m, nr = 1, nc = 1)
    G = makeG(phi)
    W = sigma2b*(c(1,rep(0,p-1))%*%t(c(1,rep(0,p-1))))
    C0 = makeC0(phi)
    x0 = rep(t(chol(C0))%*%rnorm(p,0,1),1)
  } else if(mod == "ar") {
    p = length(rho)
    F = array(0, c(1,p,nt)); F[1,1,] = 1
    V = matrix(sigma2m, nr = 1, nc = 1)
    G = makeG(rho)
    W = sigma2s*(c(1,rep(0,p-1))%*%t(c(1,rep(0,p-1))))
    C0 = makeC0(rho)
    x0 = rep(t(chol(C0))%*%rnorm(p,0,1),1)
  } else if(mod == "both") {
    p1 = length(phi); p2 = length(rho)
    F = array(0, c(1,p1+p2,nt)); F[1,1,] = conv; F[1,1+p1,] = 1
    V = matrix(sigma2m, nr = 1, nc = 1)
    G = bdiag(list(makeG(phi),makeG(rho)))
    eb = c(sigma2b,rep(0,p1-1))
    es = c(sigma2s,rep(0,p2-1))
    W = bdiag(list(eb%*%t(eb),es%*%t(es)))
    C0 = bdiag(list(makeC0(phi),makeC0(rho)))
    x0 = rep(t(chol(C0))%*%rnorm(p1+p2,0,1),1)
  } else { stop("mod must be 'dr','ar', or 'both'")}

  mysims = list()
  for(j in 1:N) mysims[[j]] = dlm.sim(nt, F, G, V, W, x0, beta, U)
  
  # Return data
  return(mysims)
}

# Simulate dynamic regression models with different values of phi, sigma2b
N = 30
M101_sim <- function(nsims, phi, sigma2b) dlm_ar_sim(nsims, "dr", beta = c(900,1), sigma2m = 1, phi = phi, sigma2b = sigma2b)
mydata = expand.grid(phi = seq(.1,.9,.2), sigma2b = seq(1,9,2), n = N)
M101_dat <- list(mlply(mydata, M101_sim), mydata)
save(M101_dat, file=paste(dpath,"dlm_ar_sim-",N,"-M101.rdata",sep=""))

# Simulate regression models with AR(1) errors for different values of rho, sigma2s
M010_sim <- function(nsims, rho, sigma2s) dlm_ar_sim(nsims, "ar", beta = c(900,1), rho = rho, sigma2s = sigma2s)
mydata = expand.grid(rho = seq(.1,.9,.2), sigma2s = seq(1,9,2), n = N)
M010_dat <- list(mlply(mydata, M010_sim), mydata)
save(M010_dat, file=paste(dpath,"dlm_ar_sim-",N,"-M010.rdata",sep=""))

# Simulate regression models with AR(1)+WN errors for different values of rho, sigma2s, sigma2m
M011_sim <- function(nsims, rho, sigma2s) dlm_ar_sim(nsims, "ar", beta = c(900,1), sigma2m = 1, rho = rho, sigma2s = sigma2s)
mydata = expand.grid(rho = seq(.1,.9,.2), sigma2s = seq(1,9,2),n = N)
M011_dat <- list(mlply(mydata, M011_sim), mydata)
save(M011_dat, file=paste(dpath,"dlm_ar_sim-",N,"-M011.rdata",sep=""))

# Simulate regression models with AR(2) errors for different values of rho1, sigma2s
M020_sim <- function(nsims, rho1, rho2, sigma2s) dlm_ar_sim(nsims, "ar", beta = c(900,1), rho = c(rho1, .08), sigma2s = sigma2s)
mydata = expand.grid(rho1 = seq(.1,.9,.2), sigma2s = seq(1,9,2), n = N)
M020_dat <- list(mlply(mydata, M020_sim), mydata)
save(M020_dat, file=paste(dpath,"dlm_ar_sim-",N,"-M020.rdata",sep=""))

# Simulate dynamic regression + AR(1) error models with different values of phi, sigma2b, rho, sigma2s
M111_sim <- function(nsims, phi, sigma2b, rho, sigma2s) dlm_ar_sim(nsims, "both", beta = c(900,1), sigma2m = 1, phi = phi, sigma2b = sigma2b, rho = rho, sigma2s = sigma2s)
mydata = expand.grid(phi = seq(.5,.9,.2), sigma2b = seq(1,9,4), rho = seq(.3,.7,.2), sigma2s = seq(1,9,4), n = N)
M111_dat <- list(mlply(mydata, M111_sim), mydata)
save(M111_dat, file=paste(dpath,"dlm_ar_sim-",N,"-M111.rdata",sep=""))

# Plot data
dlm_ar_plot <- function(N, n, n.sim, label)
{
  load(paste(dpath,"dlm_ar_sim-",N,"-",label,".rdata",sep=""))
  if(label == "M101")
  {
    mysim = M101_dat[[1]][[n]][[n.sim]]
    s2m = mysim$true.params$V[1,1]
    title = expression(paste("Simulated time series from ",M[101],sep=""))
    subtext = substitute(paste(phi," = ",a1,", ",sigma[b]^2," = ",a2,", ",sigma[m]^2," = ",a3,sep=""),list(a1 = M101_dat[[2]]$phi[n], a2 = M101_dat[[2]]$sigma2b[n], a3 = s2m))
  } else if(label == "M010") {
    mysim = M010_dat[[1]][[n]][[n.sim]]
    title = expression(paste("Simulated time series from ",M['010'],sep=""))
    subtext = substitute(paste(rho," = ",a1,", ",sigma[s]^2," = ",a2,sep=""),list(a1 = M010_dat[[2]]$rho[n], a2 = M010_dat[[2]]$sigma2s[n]))
  } else if(label == "M011") {
    mysim = M011_dat[[1]][[n]][[n.sim]]
    s2m = mysim$true.params$V[1,1]
    title = expression(paste("Simulated time series from ",M['011'],sep=""))
    subtext = substitute(paste(rho," = ",a1,", ",sigma[s]^2," = ",a2,", ",sigma[m]^2," = ",a3,sep=""),list(a1 = M011_dat[[2]]$rho[n], a2 = M011_dat[[2]]$sigma2s[n], a3 = s2m))
  } else if(label == "M020") {
    mysim = M020_dat[[1]][[n]][[n.sim]]
    rho2 = mysim$true.params$G[2,1]
    title = expression(paste("Simulated time series from ",M['020'],sep=""))
    subtext = substitute(paste(rho[1]," = ",a1,", ",rho[2]," = ",a2,", ",sigma[s]^2," = ",a3,sep=""),list(a1 = M020_dat[[2]]$rho[n], a2 = rho2, a3 = M020_dat[[2]]$sigma2s[n]))
  } else if(label == "M111") {
    mysim = M111_dat[[1]][[n]][[n.sim]]
    s2m = mysim$true.params$V[1,1]
    title = expression(paste("Simulated time series from ",M[111],sep=""))
    subtext = substitute(paste(phi," = ",a1,", ",sigma[b]^2," = ",a2,", ",rho, " = ",a3,", ",sigma[s]^2," = ",a4,", ",sigma[m]^2," = ",a5,sep=""),list(a1 = M111_dat[[2]]$phi[n], a2 = M111_dat[[2]]$sigma2b[n], a3 = M111_dat[[2]]$rho[n], a4 = M111_dat[[2]]$sigma2s[n], a5 = s2m))
  }

  nt = dim(mysim$y)[2]
  pdf(file=paste(gpath,"dlm_ar_sim-",N,"-",label,"-",n,"-",n.sim,".pdf",sep=""))
  par(mfrow=c(3,1),mar=c(3,5,4,2)+0.1)
  plot(0:nt,c(NA,mysim$y[1,]),type="l",main=title,xlab="",ylab=expression(y[t]),cex.lab=1.5,cex.main=1.5)
  mtext(subtext, side = 3, cex = 0.85)
  par(mar = c(5,5,2,2)+0.1)
  plot(0:nt,c(NA,mysim$true.params$U[1,2,]),type="l",xlab="",ylab=expression(conv[t]),cex.lab=1.5)
  par(mar = c(7,5,0,2)+0.1)
  if(label == "M111"){
    ymin = min(mysim$x); ymax = max(mysim$x)
    plot(0:nt,mysim$x[1,],type="l",col=2,ylim = c(ymin,ymax),xlab=expression(t),ylab=expression(x[t]),cex.lab=1.5)     
    lines(0:nt,mysim$x[2,],col=4)
    legend("topleft",c(expression(beta[0]),expression(beta[1])),lty=c(1,1),col=c(4,2))
  } else {
    plot(0:nt,mysim$x[1,],type="l",xlab=expression(t),ylab=expression(x[t]),cex.lab=1.5)
  }
  dev.off()
}

require(plyr)
#mydata = data.frame(N=rep(30,5),n=rep(5,5),n.sim=rep(1,5),label=c("M101","M010","M011","M020","M111"),stringsAsFactors=FALSE)
data1 = expand.grid(N = 30, n = 1:25, n.sim = 1, label = c("M101","M010","M011","M020"),stringsAsFactors=FALSE)
data2 = expand.grid(N = 30, n = 1:81, n.sim = 1, label = "M111", stringsAsFactors=FALSE)
mydata = rbind(data1,data2)
m_ply(mydata, dlm_ar_plot)
