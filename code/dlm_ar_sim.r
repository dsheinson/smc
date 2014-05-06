source("dlm_sim.r")
source("dlm_ar_functions.r")

# Set data and graphics path
dpath = "../data/"
gpath = "../graphs/"

# Function to simulate data from dlm for different regression models
# For dynamic regression models, assume first 'nd' non-intercept regression coefficients are dynamic (rest are constant)
dlm_ar_sim <- function(N, mod, beta, nd = 1, sigma2m = 0, phi = NULL, sigma2b = 0, rho = NULL, sigma2s = 0)
{
  beta = rep(beta,1); stopifnot(length(beta) == dim(X)[2])
  d = length(beta); stopifnot(d > 1)
  U = array(t(X), c(1,d,nt))

  if(mod == "dr") 
  {
    stopifnot(is.numeric(phi) & length(sigma2b) == 1 & nd < d)
    p = length(phi)
    Xd = U[1,2:(1+nd),]
    if(!is.matrix(Xd)) Xd = t(Xd)
    F = array(0, c(1,p,nt)); F[1,1,] = apply(Xd, 2, sum)
    V = matrix(sigma2m, nr = 1, nc = 1)
    G = makeG(phi)
    W = sigma2b*(c(1,rep(0,p-1))%*%t(c(1,rep(0,p-1))))
    C0 = makeC0(phi)
    x0 = rep(t(chol(C0))%*%rnorm(p,0,1),1)
  } else if(mod == "dr2") {
    stopifnot(is.list(phi) & length(sigma2b) > 1 & length(sigma2b) < d)
    stopifnot(length(phi) == length(sigma2b))
    p = rep(NA, length(phi))
    for(i in 1:length(phi)) p[i] = length(phi[[i]])
    F = array(0, c(1,sum(p),nt)); F[1,1,] = U[1,2,]
    for(i in 2:length(p)) F[1,1+sum(p[1:(i-1)]),] = U[1,i+1,] 
    V = matrix(sigma2m, nr = 1, nc = 1)
    G = makeG(phi[[1]])
    for(i in 2:length(p)) G = bdiag(list(G,makeG(phi[[i]])))
    W = matrix(0,nr=sum(p),nc=sum(p)); W[1,1] = sigma2b[1]
    for(i in 2:length(p)) W[1+sum(p[1:(i-1)]),1+sum(p[1:(i-1)])] = sigma2b[i]
    C0 = makeC0(phi[[1]])
    for(i in 2:length(p)) C0 = bdiag(list(C0,makeC0(phi[[i]])))
    x0 = rep(t(chol(C0))%*%rnorm(sum(p),0,1),1)
  } else if(mod == "ar") {
    p = length(rho)
    F = array(0, c(1,p,nt)); F[1,1,] = 1
    V = matrix(sigma2m, nr = 1, nc = 1)
    G = makeG(rho)
    W = sigma2s*(c(1,rep(0,p-1))%*%t(c(1,rep(0,p-1))))
    C0 = makeC0(rho)
    x0 = rep(t(chol(C0))%*%rnorm(p,0,1),1)
  } else if(mod == "both") {
    stopifnot(is.numeric(phi) & length(sigma2b) == 1 & nd < d)
    p1 = length(phi); p2 = length(rho)
    Xd = U[1,2:(2+nd),]
    if(!is.matrix(Xd)) Xd = t(Xd)
    F = array(0, c(1,p1+p2,nt)); F[1,1,] = apply(Xd, 2, sum); F[1,1+p1,] = 1
    V = matrix(sigma2m, nr = 1, nc = 1)
    G = bdiag(list(makeG(phi),makeG(rho)))
    eb = c(sigma2b,rep(0,p1-1))
    es = c(sigma2s,rep(0,p2-1))
    W = bdiag(list(eb%*%t(eb),es%*%t(es)))
    C0 = bdiag(list(makeC0(phi),makeC0(rho)))
    x0 = rep(t(chol(C0))%*%rnorm(p1+p2,0,1),1)
  } else if(mod == "both2") {
    stopifnot(is.list(phi) & length(sigma2b) > 1 & length(sigma2b) < d)
    stopifnot(length(phi) == length(sigma2b))
    p1 = rep(NA, length(phi))
    for(i in 1:length(phi)) p1[i] = length(phi[[i]])
    p2 = length(rho)
    F = array(0, c(1,sum(p1)+p2,nt)); F[1,1,] = U[1,2,]
    for(i in 2:length(p1)) F[1,1+sum(p1[1:(i-1)]),] = U[1,i+1,] 
    F[1,1+sum(p1),] = 1
    V = matrix(sigma2m, nr = 1, nc = 1)
    G = makeG(phi[[1]])
    for(i in 2:length(p1)) G = bdiag(list(G,makeG(phi[[i]])))
    G = bdiag(list(G,makeG(rho)))    
    W = matrix(0,nr=sum(p1)+p2,nc=sum(p1)+p2); W[1,1] = sigma2b[1]
    for(i in 2:length(p1)) W[1+sum(p1[1:(i-1)]),1+sum(p1[1:(i-1)])] = sigma2b[i]
    W[1+sum(p1),1+sum(p1)] = sigma2s
    C0 = makeC0(phi[[1]])
    for(i in 2:length(p1)) C0 = bdiag(list(C0,makeC0(phi[[i]])))
    C0 = bdiag(list(C0,makeC0(rho)))
    x0 = rep(t(chol(C0))%*%rnorm(sum(p1)+p2,0,1),1)
  } else { stop("mod must be 'dr', 'dr2', 'ar', 'both', or 'both2'")}

  mysims = list()
  for(j in 1:N) mysims[[j]] = dlm.sim(nt, F, G, V, W, x0, beta, U)
  
  # Return data
  cat(mod,beta,"\n")
  return(mysims)
}

# Load design matrix for single event experiment
load(paste(dpath,"fmri-design-1.rdata",sep=""))
X = fmri.design$X
nt = dim(X)[1] # assume X includes first column for intercept

set.seed(45)

# Simulate dynamic regression models with different values of phi, sigma2b
require(plyr)
N = 20
M101_sim <- function(nsims, phi, sigma2b) dlm_ar_sim(nsims, "dr", beta = c(900,5), nd = 1, sigma2m = 1, phi = phi, sigma2b = sigma2b)
mydata = expand.grid(nsims = N, phi = c(.7,.8,.9,.99), sigma2b = seq(.5,2,.5))
M101_dat <- list(mlply(mydata, M101_sim), mydata)
save(M101_dat, file=paste(dpath,"dlm_ar_sim-",N,"-M101-",dim(X)[2],".rdata",sep=""))

# Simulate regression models with AR(1) errors for different values of rho, sigma2s
M010_sim <- function(nsims, rho, sigma2s) dlm_ar_sim(nsims, "ar", beta = c(900,5), rho = rho, sigma2s = sigma2s)
mydata = expand.grid(nsims = N, rho = c(.7,.8,.9,.99), sigma2s = seq(.5,2,.5))
M010_dat <- list(mlply(mydata, M010_sim), mydata)
save(M010_dat, file=paste(dpath,"dlm_ar_sim-",N,"-M010-",dim(X)[2],".rdata",sep=""))

# Simulate regression models with AR(1)+WN errors for different values of rho, sigma2s, sigma2m
M011_sim <- function(nsims, rho, sigma2s) dlm_ar_sim(nsims, "ar", beta = c(900,5), sigma2m = 1, rho = rho, sigma2s = sigma2s)
mydata = expand.grid(nsims = N, rho = c(.7,.8,.9,.99), sigma2s = seq(.5,2,.5))
M011_dat <- list(mlply(mydata, M011_sim), mydata)
save(M011_dat, file=paste(dpath,"dlm_ar_sim-",N,"-M011-",dim(X)[2],".rdata",sep=""))

# Simulate regression models with AR(2) errors for different values of rho2, sigma2s
M020_sim <- function(nsims, rho2, sigma2s) dlm_ar_sim(nsims, "ar", beta = c(900,5), rho = c(.6, rho2), sigma2s = sigma2s)
mydata = expand.grid(nsims = N, rho2 = c(-.3,-.1,.1,.3), sigma2s = seq(.5,2,.5))
M020_dat <- list(mlply(mydata, M020_sim), mydata)
save(M020_dat, file=paste(dpath,"dlm_ar_sim-",N,"-M020-",dim(X)[2],".rdata",sep=""))

# # Simulate dynamic regression + AR(1) error models with different values of phi, sigma2b, rho, sigma2s
# M111_sim <- function(nsims, phi, sigma2b, rho, sigma2s) dlm_ar_sim(nsims, "both", beta = c(900,5), nd = 1, sigma2m = .5, phi = phi, sigma2b = sigma2b, rho = rho, sigma2s = sigma2s)
# mydata = expand.grid(nsims = N, phi = c(.9,.95,.99), sigma2b = seq(.5,2,.5), rho = seq(.3,.7,.2), sigma2s = seq(.5,2,.5))
# M111_dat <- list(mlply(mydata, M111_sim), mydata)
# save(M111_dat, file=paste(dpath,"dlm_ar_sim-",N,"-M111-",dim(X)[2],".rdata",sep=""))

# Load design matrix for 2-event experiment
load(paste(dpath,"fmri-design-2.rdata",sep=""))
X = fmri.design$X
nt = dim(X)[1] # assume X includes first column for intercept

# Set seed
set.seed(61)

# Simulate dynamic regression models (same beta_t's) with different values of phi, sigma2b
require(plyr)
N = 20
M101_sim <- function(nsims, phi, sigma2b) dlm_ar_sim(nsims, "dr", beta = c(900,5,2), nd = 2, sigma2m = 1, phi = phi, sigma2b = sigma2b)
mydata = expand.grid(nsims = N, phi = c(.7,.9,.95,.99), sigma2b = seq(.5,2,.5))
M101_dat <- list(mlply(mydata, M101_sim), mydata)
save(M101_dat, file=paste(dpath,"dlm_ar_sim-",N,"-M101-",dim(X)[2],".rdata",sep=""))

# Simulate dynamic regression models (diff beta_t's) with different values of phi1, phi2, sigma2b
M101_sim <- function(nsims, phi1, phi2, sigma2b1, sigma2b2)
{
  phi = list(phi1,phi2); sigma2b = c(sigma2b1,sigma2b2)
  dlm_ar_sim(nsims, "dr2", beta = c(900,5,2), nd = 2, sigma2m = 1, phi = phi, sigma2b = sigma2b)
}
mydata = expand.grid(nsims = N, phi1 = c(.7,.99), phi2 = c(.7,.99), sigma2b1 = c(.5,1.5), sigma2b2 = c(.5,1.5))
M101_dat <- list(mlply(mydata, M101_sim), mydata)
save(M101_dat, file=paste(dpath,"dlm_ar_sim-",N,"-M101-diff-",dim(X)[2],".rdata",sep=""))

# Simulate regression models with AR(1) errors for different values of rho, sigma2s
M010_sim <- function(nsims, rho, sigma2s) dlm_ar_sim(nsims, "ar", beta = c(900,5,2), rho = rho, sigma2s = sigma2s)
mydata = expand.grid(rho = seq(.2,.8,.2), sigma2s = seq(.5,2,.5), n = N)
M010_dat <- list(mlply(mydata, M010_sim), mydata)
save(M010_dat, file=paste(dpath,"dlm_ar_sim-",N,"-M010-",dim(X)[2],".rdata",sep=""))

# Simulate regression models with AR(1)+WN errors for different values of rho, sigma2s, sigma2m
M011_sim <- function(nsims, rho, sigma2s) dlm_ar_sim(nsims, "ar", beta = c(900,5,2), sigma2m = .5, rho = rho, sigma2s = sigma2s)
mydata = expand.grid(rho = seq(.2,.8,.2), sigma2s = seq(.5,2,.5),n = N)
M011_dat <- list(mlply(mydata, M011_sim), mydata)
save(M011_dat, file=paste(dpath,"dlm_ar_sim-",N,"-M011-",dim(X)[2],".rdata",sep=""))

# Simulate regression models with AR(2) errors for different values of rho1, sigma2s
M020_sim <- function(nsims, rho1, rho2, sigma2s) dlm_ar_sim(nsims, "ar", beta = c(900,5,2), rho = c(rho1, .08), sigma2s = sigma2s)
mydata = expand.grid(rho1 = seq(.2,.8,.2), sigma2s = seq(.5,2,.5), n = N)
M020_dat <- list(mlply(mydata, M020_sim), mydata)
save(M020_dat, file=paste(dpath,"dlm_ar_sim-",N,"-M020-",dim(X)[2],".rdata",sep=""))

# # Simulate dynamic regression + AR(1) error models (same beta_t's) with different values of phi, sigma2b, rho, sigma2s
# M111_sim <- function(nsims, phi, sigma2b, rho, sigma2s) dlm_ar_sim(nsims, "both", beta = c(900,5,2), nd = 2, sigma2m = .5, phi = phi, sigma2b = sigma2b, rho = rho, sigma2s = sigma2s)
# mydata = expand.grid(phi = c(.9,.95,.99), sigma2b = seq(.5,2,.5), rho = seq(.3,.7,.2), sigma2s = seq(.5,2,.5), n = N)
# M111_dat <- list(mlply(mydata, M111_sim), mydata)
# save(M111_dat, file=paste(dpath,"dlm_ar_sim-",N,"-M111-",dim(X)[2],".rdata",sep=""))

# Plot data
dlm_ar_plot <- function(N, n, n.sim, label, d)
{
  load(paste(dpath,"dlm_ar_sim-",N,"-",label,"-",d,".rdata",sep=""))
  if(label == "M101")
  {
    mysim = M101_dat[[1]][[n]][[n.sim]]
    s2m = mysim$true.params$V[1,1]
    title = expression(paste("Simulated time series from ",M[101],sep=""))
    subtext = substitute(paste(beta[0]," = ",b0,", ",beta[1]," = ",b1,", ",beta[2]," = ",b2,", ",phi," = ",a1,", ",sigma[b]^2," = ",a2,", ",sigma[m]^2," = ",a3,sep=""),list(b0 = mysim[[3]]$beta[1], b1 = mysim[[3]]$beta[2], b2 = mysim[[3]]$beta[3], a1 = M101_dat[[2]]$phi[n], a2 = M101_dat[[2]]$sigma2b[n], a3 = s2m))
  } else if(label == "M101-diff") {
    mysim = M101_dat[[1]][[n]][[n.sim]]
    s2m = mysim$true.params$V[1,1]
    title = expression(paste("Simulated time series from ",M[101],sep=""))
    subtext = substitute(paste(beta[0]," = ",b0,", ",beta[1]," = ",b1,", ",beta[2]," = ",b2,", ",phi['1,1']," = ",a1,", ",phi['2,1']," = ",a2,", ",sigma['b,1']^2," = ",a3,", ",sigma['b,2']^2," = ",a4,", ",sigma[m]^2," = ",a5,sep=""),list(b0 = mysim[[3]]$beta[1], b1 = mysim[[3]]$beta[2], b2 = mysim[[3]]$beta[3], a1 = M101_dat[[2]]$phi1[n], a2 = M101_dat[[2]]$phi2[n], a3 = M101_dat[[2]]$sigma2b1[n], a4 = M101_dat[[2]]$sigma2b2[n], a5 = s2m))
  } else if(label == "M010") {
    mysim = M010_dat[[1]][[n]][[n.sim]]
    title = expression(paste("Simulated time series from ",M['010'],sep=""))
    subtext = substitute(paste(beta[0]," = ",b0,", ",beta[1]," = ",b1,", ",beta[2]," = ",b2,", ",rho," = ",a1,", ",sigma[s]^2," = ",a2,sep=""),list(b0 = mysim[[3]]$beta[1], b1 = mysim[[3]]$beta[2], b2 = mysim[[3]]$beta[3], a1 = M010_dat[[2]]$rho[n], a2 = M010_dat[[2]]$sigma2s[n]))
  } else if(label == "M011") {
    mysim = M011_dat[[1]][[n]][[n.sim]]
    s2m = mysim$true.params$V[1,1]
    title = expression(paste("Simulated time series from ",M['011'],sep=""))
    subtext = substitute(paste(beta[0]," = ",b0,", ",beta[1]," = ",b1,", ",beta[2]," = ",b2,", ",rho," = ",a1,", ",sigma[s]^2," = ",a2,", ",sigma[m]^2," = ",a3,sep=""),list(b0 = mysim[[3]]$beta[1], b1 = mysim[[3]]$beta[2], b2 = mysim[[3]]$beta[3], a1 = M011_dat[[2]]$rho[n], a2 = M011_dat[[2]]$sigma2s[n], a3 = s2m))
  } else if(label == "M020") {
    mysim = M020_dat[[1]][[n]][[n.sim]]
    rho2 = mysim$true.params$G[2,1]
    title = expression(paste("Simulated time series from ",M['020'],sep=""))
    subtext = substitute(paste(beta[0]," = ",b0,", ",beta[1]," = ",b1,", ",beta[2]," = ",b2,", ",rho[1]," = ",a1,", ",rho[2]," = ",a2,", ",sigma[s]^2," = ",a3,sep=""),list(b0 = mysim[[3]]$beta[1], b1 = mysim[[3]]$beta[2], b2 = mysim[[3]]$beta[3], a1 = M020_dat[[2]]$rho[n], a2 = rho2, a3 = M020_dat[[2]]$sigma2s[n]))
   } #else if(label == "M111") {
#     mysim = M111_dat[[1]][[n]][[n.sim]]
#     s2m = mysim$true.params$V[1,1]
#     title = expression(paste("Simulated time series from ",M[111],sep=""))
#     subtext = substitute(paste(phi," = ",a1,", ",sigma[b]^2," = ",a2,", ",rho, " = ",a3,", ",sigma[s]^2," = ",a4,", ",sigma[m]^2," = ",a5,sep=""),list(a1 = M111_dat[[2]]$phi[n], a2 = M111_dat[[2]]$sigma2b[n], a3 = M111_dat[[2]]$rho[n], a4 = M111_dat[[2]]$sigma2s[n], a5 = s2m))
#   }

  nt = dim(mysim$y)[2]
  pdf(file=paste(gpath,"dlm_ar_sim-",N,"-",label,"-",n,"-",n.sim,".pdf",sep=""))
  par(mfrow=c(3,1),mar=c(3,5,4,2)+0.1)
  plot(0:nt,c(NA,mysim$y[1,]),type="l",main=title,xlab="",ylab=expression(y[t]),cex.lab=1.5,cex.main=1.5)
  mtext(subtext, side = 3, cex = 0.85)
  par(mar = c(5,5,2,2)+0.1)
  ymax = max(mysim$true.params$U)
  plot(0:nt,c(NA,mysim$true.params$U[1,1,]),ylim=c(0,ymax),type="l",xlab="",ylab=expression(conv[t]),cex.lab=1.5)
  if(d > 1) for(j in 2:d) lines(0:nt, c(NA,mysim$true.params$U[1,j,]),col=j)
  par(mar = c(5,5,2,2)+0.1)
  if(label == "M101" | label == "M101-diff"){
    ymin = min(0,mysim$x); ymax = max(0,mysim$x)
    plot(0:nt,mysim$x[1,],type="l",col=2,ylim = c(ymin,ymax),xlab="TR",ylab=expression(x[t]),cex.lab=1.5)     
    abline(h=0)
    if(dim(mysim$x)[1] > 1)
    {
      lines(0:nt,mysim$x[2,],col=3)
      legend("topleft",expression(beta[0],beta[1],beta[2]),lty=rep(1,3),col=c(1,2,3))
    } else {legend("topleft",expression(beta[0],beta[1]),lty=rep(1,2),col=c(1,2))}
  } else {
    plot(0:nt,mysim$x[1,],type="l",xlab="TR",ylab=expression(x[t]),cex.lab=1.5)
  }
  dev.off()
}

require(plyr)
mydata = expand.grid(N = 20, n = 1:16, n.sim = 1, label = c("M101","M010","M011","M020"),d=2,stringsAsFactors=FALSE)
m_ply(mydata, dlm_ar_plot)
