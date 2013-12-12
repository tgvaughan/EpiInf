# Plot particle trajectories

library(rjson)
df <- fromJSON(file='SMCdebug.json')

# find bounds
intmax <- length(df)
tmax <- max(df[[intmax]]$p0$t)
nParticles <- length(df[[1]])

Imax <- 0
for (idx in 1:intmax) {
    for (p in 1:nParticles) {
        Imax <- max(Imax, df[[idx]][[p]]$n)
    }
}

plot(df[[1]][[1]]$t, df[[1]][[1]]$n, xlim=c(0,tmax), ylim=c(0,Imax), 's')
for (idx in 1:intmax) {
    for (p in 1:nParticles) {
        lines(df[[idx]][[p]]$t, df[[idx]][[p]]$n, 's')
    }
}
