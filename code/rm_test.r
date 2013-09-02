# Simulate data from a dlm with common state/observation variance
v = 1; F = 1; G = 1; V = 1; W = 1
x0 = rnorm(1,0,sqrt(v))
n = 200
x = rep(NA,n+1)
y = rep(NA,n)
x[1] = x0
for(i in 1:n)
{
  x[i+1] = G*x[i] + W*rnorm(1,0,v)
  y[i] = F*x[i] + V*rnorm(1,0,v)
}

# Plot data
pdf("../graphs/rm_test-data.pdf")
gmin = min(x,y); gmax = max(x,y)
plot(0:n,x,ylim=c(gmin,gmax),type="l",xlab=expression(t),ylab="Position")
points(1:n,y)
legend("bottomright",legend=expression(x,y),lty=c(1,NA),pch=c(NA,1))
dev.off()

# Compute true marginal filtered distributions of states and unknown variance
m = C = a = b = rep(NA,n+1)
m[1] = 0; C[1] = 1; a[1] = 1; b[1] = 1
for(i in 1:n)
{
  A = G*m[i]; R = G*C[i]*G + W
  f = F*A; Q = F*R*F + V  
  e = y[i] - f; Qinv = 1/Q
  RFQinv = R*F*(1/Q)
  m[i+1] = A + RFQinv*e
  C[i+1] = R - RFQinv*F*R
  a[i+1] = a[i] + 1/2
  b[i+1] = b[i] + (1/2)*e*Qinv*e
}

# Save image of simulated data and true posterior distribution
save.image("../data/rm_test-truth.rdata")
#load("../data/rm_test-truth.rdata")

# Function to run a particle filter
pf_run <- function(pf, np, a, b, label = "", ...){
  if(pf == "KD"){
    dllik <- function(y, x, theta) dnorm(y,x,sqrt(exp(theta)),log=TRUE)
    revo <- function(x, theta) rnorm(1,x,sqrt(exp(theta)))
    rprior <- function(j)
    {
      mytheta = 1 / rgamma(1,a,b)
      mystate = rnorm(1,0,sqrt(mytheta))
      return(list(x=mystate,theta=log(mytheta)))
    }
    pstate <- function(x, theta) x
    source("kd_pf.r")
    out.kd = kd_pf(y, dllik, pstate, revo, rprior, np, ...)
    save(out.kd,file=paste("../data/rm_test-KD-",np,"-",a,"-",b,"-",label,".rdata",sep=""))
  } else if(pf == "RM1") {
    source("rm_test_functions.r")
    rprior1 <- function(j) rprior(j,a,b)
    rmove <- function(y, x, theta) return(list(state=x,theta=sample.theta(y, x, theta)))
    source("rm_pf.r")
    out.rm1 = rm_pf(y, dllik, revo, rprior1, rmove, np, ...)
    save(out.rm1,file=paste("../data/rm_test-RM1-",np,"-",a,"-",b,"-",label,".rdata",sep=""))
  } else if(pf == "RM2") {
    source("rm_test_functions.r")
    rprior1 <- function(j) rprior(j,a,b)
    rmove <- function(y, x, theta) rm_mcmc(y, x, theta, 1)
    source("rm_pf.r")
    out.rm2 = rm_pf(y, dllik, revo, rprior1, rmove, np, ...)
    save(out.rm2,file=paste("../data/rm_test-RM2-",np,"-",a,"-",b,"-",label,".rdata",sep=""))
  }
}

# Run particle filters
mydata = expand.grid(pf=c("KD","RM1","RM2"), a=1, b=c(1,.25), label=seq(1,9,1), np=c(100,1000), progress = FALSE, method="stratified", nonuniformity="ess", threshold=0.8, log=FALSE, stringsAsFactors=FALSE)
require(plyr)
m_ply(.data = mydata, .fun = pf_run)