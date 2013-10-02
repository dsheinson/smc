source("mc_functions.r")
source("rm_pf.r")
source("rm_test_functions.r")

nt = 200
F <- G <- V <- W <- sigma <- 1
F.test <- c(1, 10)

N = 100
post.lprob = mod.lpost = array(NA, dim=c(N,length(F.test),2))
mod.priors = rep(1/length(F.test),length(F.test))
x = matrix(NA, nr = N, nc = nt + 1)
y = matrix(NA, nr = N, nc = nt)
for(j in 1:N)
{
  sim = dr.sim(nt, F, G, V, W, sigma)
  x[j,] = sim$x[,1]
  y[j,] = sim$y[,1]
  for(i in 1:length(F.test))
  { 
    post = dr.post(sim$y, F.test[i], G, V, W, 1, 1, 0, 1e6)
    post.lprob[j,i,1] = dr.prob(sim$y[,1], post$f[,1], post$Q[1,1,], post$a, post$b)
  }
  mod.lpost[j,,1] = post.lprob[j,,1] + log(mod.priors)  - log(sum(exp(post.lprob[j,,1])*mod.priors))
  print(j)
}

mod.post = exp(mod.lpost)
hist(mod.post[,1,1])
hist(mod.post[,2,1],add=TRUE,col=2)


#  a0 = b0 = 1
#  np = 100
#  rprior1 <- function(j) rprior(j,a0,b0)
#  rmove <- function(y, x, theta) rm_mcmc(y, x, theta, a0, b0, 1)
#  out.rm2 = rm_pf(sim$y[,1], dllik, revo, rprior1, rmove, np, method="stratified", nonuniformity="ess", threshold=0.8, log=FALSE)
#  post.prob[j,2] = dr.pf.prob(out.rm2)

#gmin = min(sim1$x, sim1$y)
#gmax = max(sim1$x, sim1$y)
#plot(0:nt,sim1$x[,1],ylim=c(gmin,gmax),type="l",xlab=expression(t),ylab="Position")
#points(1:nt,sim1$y[,1])
#legend("bottomright",legend=expression(x,y),lty=c(1,NA),pch=c(NA,1))

#lk = qt(0.025,2*post1$a)*sqrt(post1$C[1,1,]*(post1$b/post1$a)) + post1$m[,1]
#uk = qt(0.975,2*post1$a)*sqrt(post1$C[1,1,]*(post1$b/post1$a)) + post1$m[,1]
#lines(0:nt,lk,col=2)
#lines(0:nt,uk,col=2)

#lm = qt(0.025,2*post1$a[-nt])*sqrt(post1$Q[1,1,]*(post1$b[-nt]/post1$a[-nt])) + post1$f[,1]
#um = qt(0.975,2*post1$a[-nt])*sqrt(post1$Q[1,1,]*(post1$b[-nt]/post1$a[-nt])) + post1$f[,1]
#lines(1:nt,lm,col=4)
#lines(1:nt,um,col=4)