source("plotTraj.R")

#pdf('SIS_rhoSamp_trajInference.pdf', width=7, height=5)
#plotTraj("SIS_rhoSamp_Inference.traj", col=rgb(0,0,0,0.2), main="SIS")
#df <- read.table("SIS_rhoSamp_Inference_truth.traj", header=T)
#origin <- 4
#lastIdx <- which(df$eventType==2)-1
#lines(origin - df$t[1:lastIdx], df$I[1:lastIdx], col='white', lwd=4)
#lines(origin - df$t[1:lastIdx], df$I[1:lastIdx], col='red', lwd=2)
#dev.off()


#pdf('SIS_psiSamp_trajInference.pdf', width=7, height=5)
plotTraj("SIS_psiSamp_Inference.traj", col=rgb(0,0,0,0.2), main="SIS", includeFinalState=F)
df <- read.table("SIS_psiSamp_Inference_truth.traj", header=T)
origin <- 4.96998
sampIdx <- which(df$eventType==3)
lastIdx <- sampIdx[length(sampIdx)]-1
lines(origin - df$t[1:lastIdx], df$I[1:lastIdx], col='white', lwd=4)
lines(origin - df$t[1:lastIdx], df$I[1:lastIdx], col='red', lwd=2)
#dev.off()

#pdf('SIR_psiSamp_trajInference.pdf', width=7, height=5)
#plotTraj("SIR_psiSamp_Inference.traj", col=rgb(0,0,0,0.2), main="SIS")
#df <- read.table("SIR_psiSamp_Inference_truth.traj", header=T)
#origin <- 4
#lastIdx <- which(df$eventType==2)-1
#lines(origin - df$t[1:lastIdx], df$I[1:lastIdx], col='white', lwd=4)
#lines(origin - df$t[1:lastIdx], df$I[1:lastIdx], col='red', lwd=2)
#dev.off()
