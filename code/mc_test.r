source("mc_functions.r")

F <- G <- V <- W <- sigma <- 1
sim = dr.sim(nt, F, G, V, W, sigma)
post = dr.post(sim$y, F, G, V, W.test[i], 1, 1, 0, 1e6)

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

nt = 200
W.test <- c(0.5, 0.75, 1, 1.25, 1.5)

N = 100
post.lprob = mod.lpost = matrix(NA, nr=N, nc=length(W.test))
mod.priors = rep(1/length(W.test),length(W.test))
x = matrix(NA, nr = N, nc = nt + 1)
y = matrix(NA, nr = N, nc = nt)
for(j in 1:N)
{
  sim = dr.sim(nt, F, G, V, W, sigma)
  x[j,] = sim$x[,1]
  y[j,] = sim$y[,1]
  for(i in 1:length(W.test))
  { 
    post = dr.post(sim$y, F, G, V, W.test[i], 1, 1, 0, 1e6)
    post.lprob[j,i] = dr.prob(sim$y[,1], post$f[,1], post$Q[1,1,], post$a, post$b)
  }
  mod.lpost[j,] = post.lprob[j,] + log(mod.priors)  - log(sum(exp(post.lprob[j,])*mod.priors))
  print(j)
}

mod.post = exp(mod.lpost)
table(apply(mod.post, 1, function(x) which(x == max(x))))

