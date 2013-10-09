nt = 150
np = 100
N = 100
load(paste("../data/mc_pf_test-",nt,"-",np,"-",N,".rdata",sep=""))

x1 = apply(out$mod.post[,,1],1,function(x) which(x == max(x)))
x2 = apply(out$mod.post[,,2],1,function(x) which(x == max(x)))
p.same = sum(x1 == x2) / N
se.p.same = sqrt(p.same*(1-p.same) / N)
c(p.same - 1.96*se.p.same, p.same + 1.96*se.p.same)

p1 = table(x1) / N
p2 = table(x2) / N
se.p1 = sqrt(p1*(1-p1) / N)
se.p2 = sqrt(p2*(1-p2) / N)
cbind(p1 - 1.96*se.p1, p1 + 1.96*se.p1)
cbind(p2 - 1.96*se.p2, p2 + 1.96*se.p2)

pdf(file="../graphs/mod-post.pdf",width=21,height=14)
par(mfcol=c(2,3),mar=c(9,8,8,0)+0.1,mgp=c(5,1,0))
xlab = c(expression(paste("p(",M[1],"|",y[1:t],")",sep="")),expression(paste("p(",M[2],"|",y[1:t],")",sep="")),expression(paste("p(",M[3],"|",y[1:t],")",sep="")))
main = c(expression(paste("Histogram of Posterior Model Probability of ",M[1],sep="")),expression(paste("Histogram of Posterior Model Probability of ",M[2],sep="")),expression(paste("Histogram of Posterior Model Probability of ",M[3],sep="")))
ylab1 = c("True Posterior Model Probability","","")
ylab2 = c("Estimated by Particle Filtering","","")
for(i in 1:3)
{
  hist(out$mod.post[,i,1],xlim=c(0,1),xlab="",ylab=ylab1[i],main=main[i],cex.lab=2.65,cex.main=2.5,cex.axis=2.5)
  hist(out$mod.post[,i,2],xlim=c(0,1),xlab=xlab[i],ylab=ylab2[i],main="",cex.lab=2.65,cex.main=2.5,cex.axis=2.5)
}
dev.off()