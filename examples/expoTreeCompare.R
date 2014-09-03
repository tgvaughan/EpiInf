df <- read.table('logLik.txt', header=T)
dfg <- read.table('logLik_Gabriel.txt', header=T)

plot(df$beta, df$logP-dfg$logP, 'o')
