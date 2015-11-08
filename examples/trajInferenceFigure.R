source("plotTraj.R")

#pdf('SIS_rhoSamp_trajInference.pdf', width=7, height=5)
#plotTraj("SIS_rhoSamp_Inference.traj", col=rgb(0,0,0,0.2), main="SIS", includeFinalState=F)
#df <- read.table("SIS_rhoSamp_Inference_truth.traj", header=T)
#origin <- 4
#lastIdx <- which(df$eventType==2)-1
#lines(origin - df$t[1:lastIdx], df$I[1:lastIdx], col='white', lwd=4)
#lines(origin - df$t[1:lastIdx], df$I[1:lastIdx], col='red', lwd=2)
#dev.off()


#pdf('SIS_psiSamp_trajInference.pdf', width=7, height=5)
plotTraj("SIS_psiSamp_Inference.traj", col=rgb(0,0,0,0.05), main="SIS", includeFinalState=T)
df <- read.table("SIS_psiSamp_Inference_truth.traj", header=T)
origin <- as.numeric(system2("phylostat", "-n SIS_psiSamp_Inference_truth.tree origin", stdout=T))
sampIdx <- which(df$eventType==3)
lastIdx <- sampIdx[length(sampIdx)]
lines(origin - df$t[1:lastIdx], df$I[1:lastIdx], 's', col='white', lwd=4)
lines(origin - df$t[1:lastIdx], df$I[1:lastIdx], 's', col='red', lwd=2)
#dev.off()


#pdf('SIR_psiSamp_trajInference.pdf', width=7, height=5)
#plotTraj("SIR_psiSamp_Inference.traj", col=rgb(0,0,0,0.05), main="SIR", includeFinalState=T, burnin=0.1)
#df <- read.table("SIR_psiSamp_Inference_truth.traj", header=T)
#origin <- as.numeric(system2("phylostat", "-n SIR_psiSamp_Inference_truth.tree origin", stdout=T))
#sampIdx <- which(df$eventType==3)
#lastIdx <- sampIdx[length(sampIdx)]
#lines(origin - df$t[1:lastIdx], df$I[1:lastIdx], 's', col='white', lwd=4)
#lines(origin - df$t[1:lastIdx], df$I[1:lastIdx], 's', col='red', lwd=2)

#df <- fromJSON(file="SIR_output.json")
#for (i in 1:length(df$trajectories)) {
#         lines(origin - df$trajectories[[i]]$t, df$trajectories[[i]]$I, 's', col=rgb(0,0,1, 0.02))
#}
#dev.off()


#pdf('BD_psiSamp_trajInference.pdf', width=7, height=5)
#plotTraj("BD_psiSamp_Inference.traj", col=rgb(0,0,0,0.05), main="Linear Birth-Death", includeFinalState=T, burnin=0.1)
#df <- read.table("BD_psiSamp_Inference_truth.traj", header=T)
#origin <- as.numeric(system2("phylostat", "-n BD_psiSamp_Inference_truth.tree origin", stdout=T))
#sampIdx <- which(df$eventType==3)
#lastIdx <- sampIdx[length(sampIdx)]
#lines(origin - df$t[1:lastIdx], df$I[1:lastIdx], 's', col='white', lwd=4)
#lines(origin - df$t[1:lastIdx], df$I[1:lastIdx], 's', col='red', lwd=2)

#df <- fromJSON(file="SIR_output.json")
#for (i in 1:length(df$trajectories)) {
#         lines(origin - df$trajectories[[i]]$t, df$trajectories[[i]]$I, 's', col=rgb(0,0,1, 0.02))
#}
#dev.off()

