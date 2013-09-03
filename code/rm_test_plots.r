# Plot data
load("../data/rm_test-truth.rdata")
pdf("../graphs/rm_test-data.pdf")
gmin = min(x,y); gmax = max(x,y)
plot(0:n,x,ylim=c(gmin,gmax),type="l",xlab=expression(t),ylab="Position")
points(1:n,y)
legend("bottomright",legend=expression(x,y),lty=c(1,NA),pch=c(NA,1))
dev.off()

# Function to construct plots for particle filters
pf_plots <- function(np, alpha, beta, xlab, ylab)
{
  iter = 9

  source("rm_test_functions.r")
  load("../data/rm_test-truth.rdata")

  pdf(paste("../graphs/rm_test-states-",np,"-",alpha,"-",beta,".pdf",sep=""), width = 15, height = 15)
  par(mfrow=c(3,3),mar=c(7,6,2,0)+.1,mgp=c(4,1,0))
  for(label in 1:iter)
  {
    if(label == 1)
    {
      xlab = expression(t)
      ylab = expression(x)
      leg = TRUE
    } else {
      xlab = ""
      ylab = ""
      leg = FALSE
    }
    load(paste("../data/rm_test-KD-",np,"-",alpha,"-",beta,"-",label,".rdata",sep=""))
    load(paste("../data/rm_test-RM1-",np,"-",alpha,"-",beta,"-",label,".rdata",sep=""))
    load(paste("../data/rm_test-RM2-",np,"-",alpha,"-",beta,"-",label,".rdata",sep=""))
    pf.states_plot(m, C, a, b, x, out.kd, out.rm1, out.rm2, xlab, ylab, leg)
  }
  dev.off()

  pdf(paste("../graphs/rm_test-states-zeroed-",np,"-",alpha,"-",beta,".pdf",sep=""), width = 15, height = 15)
  par(mfrow=c(3,3),mar=c(7,6,2,0)+.1,mgp=c(4,1,0))
  for(label in 1:iter)
  {
    if(label == 1)
    {
      xlab = expression(t)
      ylab = expression(x)
      leg = TRUE
    } else {
      xlab = ""
      ylab = ""
      leg = FALSE
    }
    load(paste("../data/rm_test-KD-",np,"-",alpha,"-",beta,"-",label,".rdata",sep=""))
    load(paste("../data/rm_test-RM1-",np,"-",alpha,"-",beta,"-",label,".rdata",sep=""))
    load(paste("../data/rm_test-RM2-",np,"-",alpha,"-",beta,"-",label,".rdata",sep=""))
    pf.states.zeroed_plot(m, C, a, b, x, out.kd, out.rm1, out.rm2, xlab, ylab, leg)
  }
  dev.off()

  pdf(paste("../graphs/rm_test-pvalues-",np,"-",alpha,"-",beta,".pdf",sep=""), width = 15, height = 15)
  par(mfrow=c(3,3),mar=c(7,6,2,0)+.1,mgp=c(4,1,0))
  for(label in 1:iter)
  {
    if(label == 1)
    {
      xlab = expression(t)
      ylab = "P-value"
      leg = TRUE
    } else {
      xlab = ""
      ylab = ""
      leg = FALSE
    }
    load(paste("../data/rm_test-KD-",np,"-",alpha,"-",beta,"-",label,".rdata",sep=""))
    load(paste("../data/rm_test-RM1-",np,"-",alpha,"-",beta,"-",label,".rdata",sep=""))
    load(paste("../data/rm_test-RM2-",np,"-",alpha,"-",beta,"-",label,".rdata",sep=""))
    pf.pvalues_plot(m, C, a, b, out.kd, out.rm1, out.rm2, xlab, ylab, leg)
  }
  dev.off()

  pdf(paste("../graphs/rm_test-precision-",np,"-",alpha,"-",beta,".pdf",sep=""), width = 15, height = 15)
  par(mfrow=c(3,3),mar=c(7,6,2,0)+.1,mgp=c(4,1,0))
  for(label in 1:iter)
  {
    if(label == 1)
    {
      xlab = expression(t)
      ylab = expression(phi)
      leg = TRUE
    } else {
      xlab = ""
      ylab = ""
      leg = FALSE
    }
    load(paste("../data/rm_test-KD-",np,"-",alpha,"-",beta,"-",label,".rdata",sep=""))
    load(paste("../data/rm_test-RM1-",np,"-",alpha,"-",beta,"-",label,".rdata",sep=""))
    load(paste("../data/rm_test-RM2-",np,"-",alpha,"-",beta,"-",label,".rdata",sep=""))
    pf.precision_plot(m, C, a, b, v, out.kd, out.rm1, out.rm2, xlab, ylab, leg)
  }
  dev.off()
}

# Construct plots
mydata = expand.grid(np=c(100, 1000), alpha=1, beta=c(1, .25), stringsAsFactors=FALSE)
require(plyr)
m_ply(.data = mydata, .fun = pf_plots)