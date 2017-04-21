parseTrajectoryString <- function(trajString) {
    strValues <- strsplit(strsplit(trajString, ",")[[1]], ":")

    nStates <- length(strValues)

    res <- list()
    res$t <- rep(0, nStates)
    res$S <- rep(0, nStates)
    res$I <- rep(0, nStates)
    res$R <- rep(0, nStates)
    res$leap <- rep(FALSE, nStates)

    for (i in 1:nStates) {
        res$t[i] <- as.numeric(strValues[[i]][1])
        res$S[i] <- as.numeric(strValues[[i]][2])
        res$I[i] <- as.numeric(strValues[[i]][3])
        res$R[i] <- as.numeric(strValues[[i]][4])
        res$leap[i] <- strValues[[i]][5] == "TL"
        res$incidence[i] <- as.numeric(strValues[[i]][6])
        res$Re[i] <- as.numeric(strValues[[i]][7])
    }

    return(res)
}

capitalize <- function(string) {
    return(paste(toupper(substring(string,1,1)), substring(string,2), sep=''))
}

plotTraj <- function(fileName=NA, dataFrame=NA, colidx=2, burninFrac=0.1,
                     col=rgb(1,0,0,0.3), add=FALSE,
                     leapCol=col,
                     xlim=NA, ylim=NA,
                     showMean=TRUE,
                     presentTime=NA,
                     timesAreCalendarYears=FALSE,
                     maxTrajectories=NA,
                     target='prevalence',
                     xlab='Age', ylab=capitalize(target), main='Trajectory distribution', ...) {

    if (is.na(fileName) && is.na(dataFrame)) {
        cat("Either fileName or df must be specified.\n")
        return()
    }
 
    if (target != "prevalence" && target != "incidence" && target != "Re") {
        cat("target must be one of 'prevalence', 'incidence' or 'Re'")
        return()
    }
       
    if (!is.na(fileName)) {
        cat("Loading data...")
        dataFrame <- read.table(fileName, header=T, as.is=T)
        cat("done.\n")
    }

    targetFun <- switch(target,
           prevalence = function(t) { return(t$I) },
           incidence = function(t) { return(t$incidence) },
           Re = function(t) { return(t$Re) })

    # Remove burnin
    nRecords <- dim(dataFrame)[1]
    df <- dataFrame[-(1:ceiling(nRecords*burninFrac)),]

    # Update record count
    nRecords <- dim(dataFrame)[1]

    # Parse trajectories
    cat("Parsing trajectories...")


    if (!is.na(maxTrajectories)) {
        step <- ceiling(nRecords/maxTrajectories)
    } else {
        step <- 1
    }

    traj <- list()
    
    maxOccupancy <- 0
    maxHeight <- 0
    thisRecord <- 0
    for (i in seq(1,nRecords,by=step)) {
        if (!is.na(df[[colidx]][i])) {
            thisRecord <- thisRecord + 1
            
            traj[[thisRecord]] <- parseTrajectoryString(df[[colidx]][i])
            maxOccupancy <- max(maxOccupancy, targetFun(traj[[thisRecord]]))
            maxHeight <- max(maxHeight, traj[[thisRecord]]$t)
        }
    }

    nValidRecords <- thisRecord

    cat("done.\n")

    # Compute mean prevalence
    if (showMean) {
        cat(paste("Computing mean", target, "..."))
        
        meanTargetTimes <- seq(maxHeight, 0, length.out=100)
        meanTarget <- rep(0, 100)

        for (i in 1:nValidRecords) {
            if (is.null(traj[[i]]))
                next

            tidx <- 1
            for (sidx in 1:length(meanTargetTimes)) {
                tSamp <- meanTargetTimes[sidx]
                while (tidx <= length(traj[[i]]$t) && traj[[i]]$t[tidx]>tSamp) {
                    tidx <- tidx + 1
                }

                meanTarget[sidx] <- meanTarget[sidx] + targetFun(traj[[i]])[tidx]
            }
        }
        for (sidx in 1:length(meanTarget)) {
            meanTarget[sidx] <- meanTarget[sidx]/nValidRecords
        }

        cat("done.\n")
    }

    getTime <- function(times) {
        if (is.na(presentTime))
            return(times)
        else {
            if (timesAreCalendarYears) {
                presentYear <- floor(presentTime)
                return(as.Date(365*(presentTime - times - presentYear),
                    origin=as.Date(paste(presentYear,"-01-01", sep=""))))
            } else {
                return(presentTime - times)
            }
        }
    }



    # Create plot if necessary
    if (!add) {
        if (length(ylim)==1 && is.na(ylim))
            ylim = c(0, maxOccupancy)

        if (length(xlim)==1 && is.na(xlim))
            xlim=sort(getTime(c(0, maxHeight)))


        plot(getTime(traj[[1]]$t), traj[[1]]$I,
             col=NA, xlim=xlim, ylim=ylim,
             xlab=xlab, ylab=ylab, main=main, ...)
    }

    cat("Plotting...")
    for (i in 1:nValidRecords) {
        indices <- which(traj[[i]]$t>=0)
        lines(getTime(traj[[i]]$t[indices]), targetFun(traj[[i]])[indices], col=col, ...)

        midx <- which.min(traj[[i]]$t[indices])
        mval <- traj[[i]]$t[indices][midx]
        lines(getTime(c(0, mval)), rep(targetFun(traj[[i]])[indices][midx],2), col=col, ...)
    }

    if (showMean) {
        lines(getTime(meanTargetTimes), meanTarget, lwd=4, col='white')
        lines(getTime(meanTargetTimes), meanTarget, lwd=2, col='black')
    }
    cat("done.\n")
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


