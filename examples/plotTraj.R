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

plotTraj <- function(fileName, header, unconditioned=NA,
                     burninFrac=0.1,
                     col=rgb(1,0,0,0.3), colUnconditioned=rgb(0,0,1,0.3),
                     xlab='Time', ylab='Prevalence', main='Trajectory distribution', ...) {
    df <- read.table(fileName, header=T, as.is=T)

    # Determine which columns contain traj data
    colidx <- which(colnames(df)==header)
    if (!is.na(unconditioned))
        colidxUC <- which(colnames(df)==unconditioned)

    # Remove burnin
    nRecords <- dim(df)[1]
    df <- df[-(1:ceiling(nRecords*burninFrac)),]

    # Update record count
    nRecords <- dim(df)[1]

    # Parse trajectories
    maxOccupancy <- 0
    maxHeight <- 0
    traj <- list()
    trajUC <- list()
    for (i in 1:nRecords) {
        if (!is.na(df[[colidx]][i])) {
            traj[[i]] <- parseTrajectoryString(df[[colidx]][i])
            maxOccupancy <- max(maxOccupancy, traj[[i]]$S)
            maxOccupancy <- max(maxOccupancy, traj[[i]]$I)
            maxOccupancy <- max(maxOccupancy, traj[[i]]$R)
            maxHeight <- max(maxHeight, traj[[i]]$t)
        }
        if (!is.na(unconditioned) && !is.na(df[[colidxUC]][i])) {
            trajUC[[i]] <- parseTrajectoryString(df[[colidxUC]][i])
        }
    }
    plot(traj[[1]]$t, traj[[1]]$I, col=NA, xlim=c(0, maxHeight), ylim=c(0, maxOccupancy), xlab=xlab, ylab=ylab, main=main, ...)

    if (length(trajUC) > 0) {
        for (i in 1:length(trajUC)) {
            indices <- which(trajUC[[i]]$t>=0)
            lines(trajUC[[i]]$t[indices], trajUC[[i]]$I[indices], 's', col=colUnconditioned, ...)

            midx <- which.min(trajUC[[i]]$t[indices])
            mval <- trajUC[[i]]$t[indices][midx]
            lines(c(0, mval), rep(trajUC[[i]]$I[indices][midx],2), col=colUnconditioned, ...)
        }
    }

    for (i in 1:length(traj)) {
        indices <- which(traj[[i]]$t>=0)
        lines(traj[[i]]$t[indices], traj[[i]]$I[indices], 's', col=col, ...)

        midx <- which.min(traj[[i]]$t[indices])
        mval <- traj[[i]]$t[indices][midx]
        lines(c(0, mval), rep(traj[[i]]$I[indices][midx],2), col=col, ...)
    }


}
