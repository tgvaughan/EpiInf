df <- read.table('conditioning.txt', header=T)

pdf('conditioningPlot.pdf', width=8, height=4)
par(mfcol=c(1,2))

plot(df$f, df$discard, 'o', log='x',
     xlab='Fraction of epidemic sampled',
     ylab='Discarded fraction',
     main='Fraction of trajectories discarded')

plot(df$f, df$ratio, 'o', log='x',
     xlab='Fraction of epidemic sampled',
     ylab=expression(T/n),
     main='Ratio of tree length to tip count')

dev.off()
