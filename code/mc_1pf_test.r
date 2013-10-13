source("mc_functions.r")
source("rm_pf.r")
source("rm_test_functions.r")

F <- G <- V <- W <- sigma <- 1
nt = 200
sim = dr.sim(nt, F, G, V, W, sigma)
post = dr.post(sim$y, F, G, V, W, 1, 1, 0, 1e6)

burn = 1
lk = qt(0.025,2*post$a)*sqrt(post$C[1,1,]*(post$b/post$a)) + post$m[,1]
uk = qt(0.975,2*post$a)*sqrt(post$C[1,1,]*(post$b/post$a)) + post$m[,1]
lm = qt(0.025,2*post$a[-nt])*sqrt(post$Q[1,1,]*(post$b[-nt]/post$a[-nt])) + post$f[,1]
um = qt(0.975,2*post$a[-nt])*sqrt(post$Q[1,1,]*(post$b[-nt]/post$a[-nt])) + post$f[,1]
gmin = min(sim$x, sim$y, lk[-(1:burn)], lm[-(1:burn)])
gmax = max(sim$x, sim$y, uk[-(1:burn)], um[-(1:burn)])

plot(0:nt,sim$x[,1],ylim=c(gmin,gmax),type="l",xlab=expression(t),ylab="Position")
points(1:nt,sim$y[,1])
lines(0:nt,lk,col=2)
lines(0:nt,uk,col=2)
lines(1:nt,lm,col=4)
lines(1:nt,um,col=4)
legend("bottomright",legend=c(expression(x,y),"95% CI","95% PI"),lty=c(1,NA,1,1),pch=c(NA,1,NA,NA),col=c(1,1,2,4))

true.marg = dr.prob(sim$y[,1], post$f[,1], post$Q[1,1,], post$a, post$b)

pf_lprob <- function(np)
{
  a0 = b0 = 1
  rprior1 <- function(j) rprior(j,a0,b0)
  rmove <- function(y, x, theta) rm_mcmc(y, x, theta, a0, b0, 1)
  out.rm2 = rm_pf(sim$y[,1], dllik, revo, rprior1, rmove, np, progress = FALSE, method="stratified", nonuniformity="ess", threshold=0.8, log=FALSE)
  print(np)
  return(dr.pf.prob(out.rm2))
}

mydata = data.frame(np = rep(c(100, 500, 1000), rep(20, 3)))
require(plyr)
myprobs = mdply(.data = mydata, .fun = pf_lprob)
save(myprobs, file = "../data/pf_margliks.rdata")