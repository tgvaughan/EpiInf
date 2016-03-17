source("plotTraj.R")
require(rjson)

#png("SIS_psiSamp_Inference.png", width=640, height=480)
#df <- read.table("SIS_psiSamp_Inference.log", header=T)
#df <- df[-(1:dim(df)[1]*.1),]
#plot(1, type='n', xlim=c(0, max(df$origin)), ylim=c(0,400),
#     xlab="Time", ylab="Prevalence", main="SIS")
#plotTraj("SIS_psiSamp_Inference.traj", "trajUnconditioned",
#         col=rgb(0,0,1,0.1), add=T)
#plotTraj("SIS_psiSamp_Inference.traj", "trajConditioned",
#         col=rgb(1,0,0,0.1), add=T)
#origin <- as.numeric(system2("phylostat", "-n SIS_psiSamp_Inference_truth.tree origin", stdout=T))
#df <- read.table("SIS_psiSamp_Inference_truth.traj", header=T)
#sampIdx <- which(df$eventType=="PSI_SAMPLE_REMOVE")
#lastIdx <- sampIdx[length(sampIdx)]
#lines(origin - df$t[1:lastIdx], df$I[1:lastIdx], 's', col='white', lwd=4)
#lines(origin - df$t[1:lastIdx], df$I[1:lastIdx], 's', col='black', lwd=2)
#
#legend('topright', inset=0.05, c("Truth", "Conditioned", "Unconditioned"), lty=1, lwd=4, col=c("black", "red", "blue"))
#dev.off()

#png("SIR_psiSamp_Inference.png", width=640, height=480)
pdf("SIR_psiSampNoSeq_Inference.pdf", width=7, height=4)
df <- read.table("SIR_psiSampNoSeq_Inference.log", header=T)
df <- df[-(1:dim(df)[1]*0.1),]
plot(1, type="n", xlim=c(0, max(df$origin)), ylim=c(0, 200),
     xlab="Time", ylab="Prevalence", main="")

plotTraj("SIR_psiSampNoSeq_Inference.traj", "trajUnconditioned",
         col=rgb(0,0,1,0.05), add=T)
plotTraj("SIR_psiSampNoSeq_Inference.traj", "trajConditioned",
         col=rgb(1,0,0,0.05), add=T)

origin <- as.numeric(system2("phylostat", "-n SIR_psiSampNoSeq_Inference_truth.tree origin", stdout=T))
df <- read.table("SIR_psiSampNoSeq_Inference_truth.traj", header=T)
sampIdx <- which(df$eventType=="PSI_SAMPLE_REMOVE")
lastIdx <- sampIdx[length(sampIdx)]
lines(origin - df$t[1:lastIdx], df$I[1:lastIdx], col='white', lwd=4)
lines(origin - df$t[1:lastIdx], df$I[1:lastIdx], col='black', lwd=2)

legend('topright', inset=0.05,
       c("Truth", "Conditioned", "Unconditioned"),
       lty=1, lwd=4,
       col=c("black", "red", "blue"),
       )
dev.off()


#png("BD_psiSamp_Inference.png", width=640, height=480)
#df <- read.table("BD_psiSamp_Inference.log", header=T)
#df <- df[-(1:dim(df)[1]*0.1),]
#plot(1, type="n", xlim=c(0, 10), ylim=c(1, 10000), log='y',
#     xlab="Time", ylab="Prevalence", main="Linear birth-death")
#
#for (i in 1:length(df$origin)) {
#    traj <- simBD(df$origin[i], df$birthRate[i], df$removalRate[i])
#    lines(traj$t, traj$I, 's', col=rgb(0,0,1,0.2))
#}
#
#plotTraj("BD_psiSamp_Inference.traj", "trajConditioned",
#         col=rgb(1,0,0,0.2), add=T)
#origin <- as.numeric(system2("phylostat", "-n BD_psiSamp_Inference_truth.tree origin", stdout=T))
#df <- read.table("BD_psiSamp_Inference_truth.traj", header=T)
#sampIdx <- which(df$eventType=="PSI_SAMPLE_REMOVE")
#lastIdx <- sampIdx[length(sampIdx)]
#lines(origin - df$t[1:lastIdx], df$I[1:lastIdx], 's', col='white', lwd=4)
#lines(origin - df$t[1:lastIdx], df$I[1:lastIdx], 's', col='black', lwd=2)
#
#legend('topright', inset=0.05, c("Truth", "Conditioned", "Unconditioned"), lty=1, lwd=4, col=c("black", "red", "blue"))
#dev.off()


#pdf('SIS_rhoSamp_trajInference.pdf', width=7, height=5)
#plotTraj("SIS_rhoSamp_Inference.traj", col=rgb(0,0,0,0.2), main="SIS", includeFinalState=F)
#df <- read.table("SIS_rhoSamp_Inference_truth.traj", header=T)
#origin <- 4
#lastIdx <- which(df$eventType==2)-1
#lines(origin - df$t[1:lastIdx], df$I[1:lastIdx], col='white', lwd=4)
#lines(origin - df$t[1:lastIdx], df$I[1:lastIdx], col='red', lwd=2)
#dev.off()

