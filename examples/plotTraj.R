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

getTime <- function(times, presentTime) {
    if (is.na(presentTime))
        return(times)
    else
        return(presentTime - times)
}

plotTraj <- function(fileName, header="trajectory", burninFrac=0.1,
                     col=rgb(1,0,0,0.3), add=FALSE,
                     xlim=NA, ylim=NA,
                     showMean=TRUE,
                     presentTime=NA,
                     xlab='Age', ylab='Prevalence', main='Trajectory distribution', ...) {
    df <- read.table(fileName, header=T, as.is=T)

    # Determine which columns contain traj data
    colidx <- which(colnames(df)==header)

    # Remove burnin
    nRecords <- dim(df)[1]
    df <- df[-(1:ceiling(nRecords*burninFrac)),]

    # Update record count
    nRecords <- dim(df)[1]

    # Parse trajectories
    maxOccupancy <- 0
    maxHeight <- 0
    traj <- list()
#    trajUC <- list()
    for (i in 1:nRecords) {
        if (!is.na(df[[colidx]][i])) {
            traj[[i]] <- parseTrajectoryString(df[[colidx]][i])
            #maxOccupancy <- max(maxOccupancy, traj[[i]]$S)
            maxOccupancy <- max(maxOccupancy, traj[[i]]$I)
            #maxOccupancy <- max(maxOccupancy, traj[[i]]$R)
            maxHeight <- max(maxHeight, traj[[i]]$t)
        }
    }

    # Compute mean prevalence
    if (showMean) {
        meanPrevTimes <- seq(maxHeight, 0, length.out=100)
        meanPrev <- rep(0, 100)

        for (i in 1:nRecords) {
            tidx <- 1
            for (sidx in 1:length(meanPrevTimes)) {
                tSamp <- meanPrevTimes[sidx]
                while (tidx <= length(traj[[i]]$t) && traj[[i]]$t[tidx]>tSamp) {
                    tidx <- tidx + 1
                }

                meanPrev[sidx] <- meanPrev[sidx] + traj[[i]]$I[tidx]
            }
        }
        for (sidx in 1:length(meanPrev)) {
            meanPrev[sidx] <- meanPrev[sidx]/nRecords
        }
    }


    # Create plot if necessary
    if (!add) {
        if (length(ylim)==1 && is.na(ylim))
            ylim = c(0, maxOccupancy)

        if (length(xlim)==1 && is.na(xlim))
            xlim=sort(getTime(c(0, maxHeight), presentTime))


        plot(getTime(traj[[1]]$t, presentTime), traj[[1]]$I,
             col=NA, xlim=xlim, ylim=ylim,
             xlab=xlab, ylab=ylab, main=main, ...)
    }

    for (i in 1:length(traj)) {
        indices <- which(traj[[i]]$t>=0)
        lines(getTime(traj[[i]]$t[indices], presentTime), traj[[i]]$I[indices], col=col, ...)

        midx <- which.min(traj[[i]]$t[indices])
        mval <- traj[[i]]$t[indices][midx]
        lines(getTime(c(0, mval), presentTime), rep(traj[[i]]$I[indices][midx],2), col=col, ...)
    }

    if (showMean) {
        lines(getTime(meanPrevTimes, presentTime), meanPrev, lwd=4, col='white')
        lines(getTime(meanPrevTimes, presentTime), meanPrev, lwd=2, col='black')
    }
}

simSIR <- function(origin, beta, gamma, N) {
    S <- N - 1
    I <- 1
    R <- 0

    tidx <- 1
    t <- origin

    while (TRUE) {

        a_infect <- beta*S[tidx]*I[tidx]
        a_recov <- gamma*I[tidx]
        a_tot <- a_infect + a_recov

        if (a_tot > 0)
            t[tidx+1] <- t[tidx] - rexp(1, a_tot)
        else
            t[tidx+1] <- -Inf

        if (t[tidx+1]<0)
            break

        if (runif(1, min=0, max=a_tot) < a_infect) {
            # Infection
            S[tidx+1] = S[tidx] - 1
            I[tidx+1] = I[tidx] + 1
            R[tidx+1] = R[tidx]
        } else {
            # Recovery
            S[tidx+1] = S[tidx]
            I[tidx+1] = I[tidx] - 1
            R[tidx+1] = R[tidx] + 1
        }

        tidx <- tidx + 1
    }

    t[tidx+1] = 0
    S[tidx+1] = S[tidx]
    I[tidx+1] = I[tidx]
    R[tidx+1] = R[tidx]

    return(data.frame(t=t, S=S, I=I, R=R))
}

simSIS <- function(origin, beta, gamma, N) {
    S <- N - 1
    I <- 1

    tidx <- 1
    t <- origin

    while (TRUE) {

        a_infect <- beta*S[tidx]*I[tidx]
        a_recov <- gamma*I[tidx]
        a_tot <- a_infect + a_recov

        if (a_tot > 0)
            t[tidx+1] <- t[tidx] - rexp(1, a_tot)
        else
            t[tidx+1] <- -Inf

        if (t[tidx+1]<0)
            break

        if (runif(1, min=0, max=a_tot) < a_infect) {
            # Infection
            S[tidx+1] = S[tidx] - 1
            I[tidx+1] = I[tidx] + 1
        } else {
            # Recovery
            S[tidx+1] = S[tidx] + 1
            I[tidx+1] = I[tidx] - 1
        }

        tidx <- tidx + 1
    }

    t[tidx+1] = 0
    S[tidx+1] = S[tidx]
    I[tidx+1] = I[tidx]

    return(data.frame(t=t, S=S, I=I))
}

simBD <- function(origin, lambda, mu) {
    I <- 1

    tidx <- 1
    t <- origin

    while (TRUE) {

        a_birth <- lambda*I[tidx]
        a_death <- mu*I[tidx]
        a_tot <- a_birth + a_death

        if (a_tot > 0)
            t[tidx+1] <- t[tidx] - rexp(1, a_tot)
        else
            t[tidx+1] <- -Inf

        if (t[tidx+1]<0)
            break

        if (runif(1, min=0, max=a_tot) < a_birth) {
            # Birth
            I[tidx+1] = I[tidx] + 1
        } else {
            # Death
            I[tidx+1] = I[tidx] - 1
        }

        tidx <- tidx + 1
    }

    t[tidx+1] = 0
    I[tidx+1] = I[tidx]

    return(data.frame(t=t, I=I))
}


