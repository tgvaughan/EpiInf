source("plotTraj.R")

plotTraj("SIS_rhoSamp_Inference.traj", col=rgb(0,0,0,0.2), main="SIS")
df <- read.table("SIS_rhoSamp_Inference_truth.traj", header=T)
origin <- 4
lastIdx <- which(df$eventType==2)-1
lines(origin - df$t[1:lastIdx], df$I[1:lastIdx], col='white', lwd=4)
lines(origin - df$t[1:lastIdx], df$I[1:lastIdx], col='red', lwd=2)
