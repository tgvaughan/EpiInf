parseTrajectoryString <- function(trajString) {
    strValues <- strsplit(strsplit(trajString, ",")[[1]], ":")

    nStates <- length(strValues)

    res <- list()
    res$t <- rep(0, nStates)
    res$S <- rep(0, nStates)
    res$I <- rep(0, nStates)
    res$R <- rep(0, nStates)

    for (i in 1:nStates) {
        res$t[i] <- as.numeric(strValues[[i]][1])
        res$S[i] <- as.numeric(strValues[[i]][2])
        res$I[i] <- as.numeric(strValues[[i]][3])
        res$R[i] <- as.numeric(strValues[[i]][4])
    }

    return(res)
}

plotTraj <- function(fileName, burninFrac=0.1, includeFinalState=TRUE, col=rgb(0,0,0,0.3), xlab='Time', ylab='Prevalence', main='Trajectory distribution', ...) {
    df <- read.table(fileName, header=T, as.is=T)

    # Remove burnin
    nRecords <- dim(df)[1]
    df <- df[-(1:ceiling(nRecords*burninFrac)),]

    # Update record count
    nRecords <- dim(df)[1]

    # Parse trajectories
    maxOccupancy <- 0
    maxHeight <- 0
    traj <- list()
    for (i in 1:nRecords) {
        if (is.na(df$traj[i]))
            next

        traj[[i]] <- parseTrajectoryString(df$traj[i])
        maxOccupancy <- max(maxOccupancy, traj[[i]]$S)
        maxOccupancy <- max(maxOccupancy, traj[[i]]$I)
        maxOccupancy <- max(maxOccupancy, traj[[i]]$R)
        maxHeight <- max(maxHeight, traj[[i]]$t)
    }
    plot(traj[[1]]$t, traj[[1]]$I, col=NA, xlim=c(0, maxHeight), ylim=c(0, maxOccupancy), xlab=xlab, ylab=ylab, main=main, ...)

    for (i in 1:nRecords) {
        if (includeFinalState) {
            lines(traj[[i]]$t, traj[[i]]$I, 's', col=col, ...)
        } else {
            n <- length(traj[[i]]$t)
            lines(traj[[i]]$t[-n], traj[[i]]$I[-n], 's', col=col, ...)
        }
    }
}
