source("mc_functions.r")
source("rm_pf.r")
source("rm_test_functions.r")

mc_pfs <- function(nt, np, N)
{
  F <- G <- V <- W <- sigma <- 1
  W.test <- c(0.5, 1, 1.5)
  revo.test <- function(x, theta, W) rnorm(1,x,sqrt(theta*W))

  post.lprob = mod.lpost = array(NA, dim=c(N,length(W.test),2))
  mod.priors = rep(1/length(W.test),length(W.test))
  x = matrix(NA, nr = N, nc = nt + 1)
  y = matrix(NA, nr = N, nc = nt)
  a0 = b0 = 1
  rprior1 <- function(j) rprior(j,a0,b0)
  rmove <- function(y, x, theta) rm_mcmc(y, x, theta, a0, b0, 1)
  for(j in 1:N)
  {
    sim = dr.sim(nt, F, G, V, W, sigma)
    x[j,] = sim$x[,1]
    y[j,] = sim$y[,1]
    for(i in 1:length(W.test))
    { 
      post = dr.post(sim$y, F, G, V, W.test[i], 1, 1, 0, 1e6)
      post.lprob[j,i,1] = dr.prob(sim$y[,1], post$f[,1], post$Q[1,1,], post$a, post$b)
      revo1 = function(x, theta) revo.test(x, theta, W.test[i])
      out.rm2 = rm_pf(sim$y[,1], dllik, revo1, rprior1, rmove, np, method="stratified", nonuniformity="ess", threshold=0.8, log=FALSE)
      post.lprob[j,i,2] = dr.pf.prob(out.rm2)
    }
    mod.lpost[j,,1] = post.lprob[j,,1] + log(mod.priors)  - log(sum(exp(post.lprob[j,,1])*mod.priors))
    mod.lpost[j,,2] = post.lprob[j,,2] + log(mod.priors)  - log(sum(exp(post.lprob[j,,2])*mod.priors))
    print(j)
  }

  mod.post = exp(mod.lpost)
  out = list(mod.post=mod.post,post.lprob=post.lprob,x=x,y=y)
  save(out, file=paste("../data/mc_pf_test-",nt,"-",np,"-",N,".rdata",sep=""))
}

require(plyr)
mydata = expand.grid(nt=150,np=c(100,200),N=100,stringsAsFactors=FALSE)
m_ply(mydata, mc_pfs)