library(plotrix) # weighted.hist used for inference counts
library(ggplot2)
library(readr)
library(stringr)

loadData <- function(fileNames, burninFrac=0.1, nSamples=NA) {

    nFiles <- length(fileNames)
    dataFrames <- list()
    for (i in 1:nFiles) {
        dataFrames[[i]] <- read_tsv(fileNames[[i]], col_types = "ic")
    }

    # Remove burnin and downsample
    for (i in 1:nFiles) {
        nRecords <- dim(dataFrames[[i]])[1]
        dataFrames[[i]] <- dataFrames[[i]][-(1:ceiling(nRecords*burninFrac)),]

        if (!is.na(nSamples)) {
            nRecords <- dim(dataFrames[[i]])[1]
            skip <- ceiling(nRecords/nSamples)

            if (skip>0)
                dataFrames[[i]] <- dataFrames[[i]][seq(1,nRecords,by=skip),]
        }
    }

    return(dataFrames);
}


parseTrajectoryString <- function(trajString) {
  strValues <- str_split(str_split(trajString, ",")[[1]], ":", simplify = TRUE)
  
  res <- list(t =  as.numeric(strValues[,1]),
              S = as.numeric(strValues[,2]),
              I = as.numeric(strValues[,3]),
              R = as.numeric(strValues[,4]),
              leap = strValues[,5] == "TL",
              incidence = as.numeric(strValues[,6]),
              Re = as.numeric(strValues[,7]))

  return(res)
}

parseTrajectories <- function(dataFrames, colidx=2) {
    traj <- list()
    nFiles <- length(dataFrames)
    
    thisRecord <- 0
    for (i in 1:nFiles) {
        for (j in 1:dim(dataFrames[[i]])[1]) {
            if (!is.na(dataFrames[[i]][[colidx]][j])) {
                thisRecord <- thisRecord + 1

                traj[[thisRecord]] <- parseTrajectoryString(dataFrames[[i]][[colidx]][j])
                traj[[thisRecord]]$fileNum <- i
            }
        }
    }

    return(traj)
}

capitalize <- function(string) {
    words <- strsplit(string, '_')
    heads <- toupper(substring(words, 1, 1))
    tails <- tolower(substring(words, 2))

    return(paste(paste0(heads,tails), collapse=' '))
}

getInterpolatedValues <- function(traj, ages, targetFun, subSample) {
    targetValues <- matrix(NA, length(traj), 100)

    for (i in seq(1,length(traj),length.out=subSample)) {
        tidx <- 1
        for (sidx in 1:length(ages)) {
            tSamp <- ages[sidx]
            while (tidx < length(traj[[i]]$t) && traj[[i]]$t[tidx]>tSamp) {
                tidx <- tidx + 1
            }

            targetValues[i,sidx] <- targetFun(traj[[i]])[tidx]
        }
    }

    return(targetValues)
}

linesBold <- function(x, y, col='black', steps=F, widthOuter=6, widthInner=2,...) {

    if (steps)
        style <- 's'
    else
        style <- 'l'

    lines(x, y, col='white', lwd=6, style)
    lines(x, y, col=col, lwd=2, style, ...)
}

plotTraj <- function(fileNames=list(), dataFrames=NULL, traj=NULL, colidx=2, burninFrac=0.1,
                     col=rgb(1,0,0,0.3), add=FALSE,
                     leapCol=col,
                     xlim=NA, ylim=NA,
                     heightLim=NA,
                     showMean=TRUE,
                     showHPD=TRUE,
                     widthOuter=6,
                     widthInner=2,
                     presentTime=NA,
                     timesAreCalendarYears=FALSE,
                     target='prevalence',
                     subSample=NA,
                     incidencePeriod=1,
                     xlab=if (is.na(presentTime)) "Age" else "Time", ylab=capitalize(target), main='Trajectory distribution', ...) {


    if (!(tolower(target) %in% c("prevalence", "scaled_prevalence", "incidence", "re"))) {
        cat("Error: target must be one of 'prevalence', 'scaled_prevalence', 'incidence' or 'Re'")
        return()
    }

    # Load data
    if (is.null(dataFrames) && is.null(traj)) {
        if (length(fileNames) == 0) {
            cat("Error: must specify at least one input file!")
            return();
        }

        cat("Loading data...")
        dataFrames <- loadData(fileNames, burninFrac)
        cat("done.\n")
    }

    # Parse trajectories
    if (is.null(traj)) {
        cat("Parsing trajectories...")
        traj <- parseTrajectories(dataFrames, colidx)
        cat("done.\n")
    }

    if (is.na(subSample))
        subSample <- length(traj)

    # Define target function for plotting
    targetFun <- switch(target,
           prevalence = function(t) { return(t$I) },
           scaled_prevalence = function(t) { return(t$I/t$S[1]*1e5) },
           incidence = function(t) { return(t$incidence/t$S*incidencePeriod) },
           Re = function(t) { return(t$Re) })

    # Identify plot boundaries
    maxHeight <- 0
    minValue <- +Inf
    maxValue <- -Inf
    for (i in seq(1,length(traj),length.out=subSample)) {
        maxHeight <- max(maxHeight, traj[[i]]$t)
        val <- targetFun(traj[[i]])
        val <- val[which(is.finite(val))]
        maxValue <- max(maxValue, val)
        minValue <- min(minValue, val)
    }

    if (!is.na(heightLim))
        maxHeight <- heightLim

    # Compute mean prevalence
    if (showMean || showHPD) {
        cat(paste0("Computing moments ", target, "..."))
        
        targetTimes <- seq(maxHeight, 0, length.out=100)
        targetValues <- getInterpolatedValues(traj, targetTimes, targetFun, subSample)

        meanTarget <- colMeans(targetValues)
        hpdTargetLow <- apply(targetValues, 2, function(x) {return(quantile(x, probs=0.025, na.rm=TRUE))})
        hpdTargetHigh <- apply(targetValues, 2, function(x) {return(quantile(x, probs=0.975, na.rm=TRUE))})

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
        if (length(ylim)==1 && is.na(ylim)) {
            if (target == 'Re')
                ylim = c(minValue, maxValue)
            else
                ylim = c(0, maxValue)
        }

        if (length(xlim)==1 && is.na(xlim))
            xlim=sort(getTime(c(0, maxHeight)))


        plot(getTime(traj[[1]]$t), targetFun(traj[[1]]),
             col=NA, xlim=xlim, ylim=ylim,
             xlab=xlab, ylab=ylab, main=main, ...)
    }

    cat("Plotting...")
    col <- rep(col, length.out=length(traj))
    for (i in seq(1,length(traj),length.out=subSample)) {
        indices <- which(traj[[i]]$t>=0)
        lines(getTime(traj[[i]]$t[indices]), targetFun(traj[[i]])[indices], col=col[traj[[i]]$fileNum], 's', ...)

        midx <- which.min(traj[[i]]$t[indices])
        mval <- traj[[i]]$t[indices][midx]
        lines(getTime(c(0, mval)), rep(targetFun(traj[[i]])[indices][midx],2), col=col[traj[[i]]$fileNum], 's', ...)
    }

    if (showMean) {
        linesBold(getTime(targetTimes), meanTarget,
                  widthOuter=widthOuter, widthInner=widthInner)
    }

    if (showHPD) {
        linesBold(getTime(targetTimes), hpdTargetLow, lty=2,
                  widthOuter=widthOuter, widthInner=widthInner)
        linesBold(getTime(targetTimes), hpdTargetHigh, lty=2,
                  widthOuter=widthOuter, widthInner=widthInner)
    }
    cat("done.\n")
}

simSIR <- function(origin, beta, gamma, N, useTL=FALSE, nLeaps=500) {
    S <- N - 1
    I <- 1
    R <- 0

    tidx <- 1
    t <- origin

    if (useTL)
        tau <- origin/nLeaps

    while (TRUE) {

        a_infect <- beta*S[tidx]*I[tidx]
        a_recov <- gamma*I[tidx]
        a_tot <- a_infect + a_recov

        if (useTL) {

            t[tidx+1] <- t[tidx] - tau

            if (t[tidx+1] < 0)
                break

            nInfect <- rpois(1, a_infect*tau)
            nRecov <- rpois(1, a_recov*tau)

            S[tidx+1] = max(S[tidx] - nInfect, 0)
            I[tidx+1] = max(I[tidx] + nInfect - nRecov, 0)
            R[tidx+1] = max(R[tidx] + nRecov, 0)


        } else {

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

        }

        tidx <- tidx + 1
    }

    t[tidx+1] = 0
    S[tidx+1] = S[tidx]
    I[tidx+1] = I[tidx]
    R[tidx+1] = R[tidx]

    return(data.frame(t=t, S=S, I=I, R=R, Re=beta*S/gamma))
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

    return(data.frame(t=t, S=S, I=I, Re=beta*S/gamma))
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

loadParamData <- function(fileNames, burninFrac=0.1, nSamples=NA) {

    nFiles <- length(fileNames)
    dataFrames <- list()
    for (i in 1:nFiles) {
        dataFrames[[i]] <- read.table(fileNames[[i]], header=T)
    }

    # remove burnin and downsample
    for (i in 1:nFiles) {
        nRecords <- dim(dataFrames[[i]])[1]
        dataFrames[[i]] <- dataFrames[[i]][-(1:ceiling(nRecords*burninFrac)),]

        if (!is.na(nSamples)) {
            nRecords <- dim(dataFrames[[i]])[1]
            skip <- ceiling(nRecords/nSamples)

            if (skip>0)
                dataFrames[[i]] <- dataFrames[[i]][seq(1,nRecords,by=skip),]
        }
    }

    return(dataFrames);
}

simTrajectories <- function(dataFrames,
                            simulationFunction=function(record) {
                                origin=record[which(startsWith(names(record), "epiOrigin"))]
                                beta=record[which(startsWith(names(record), "infectionRate"))]
                                gamma=record[which(startsWith(names(record), "recoveryRate"))]
                                N=1+record[which(startsWith(names(record), "S0"))]
                                return(simSIR(origin, beta, gamma, N))
                            }) {
    traj <- list()
    nFiles <- length(dataFrames)
    
    thisRecord <- 0
    for (i in 1:nFiles) {
        for (j in 1:dim(dataFrames[[i]])[1]) {
            thisRecord <- thisRecord + 1

            traj[[thisRecord]] <- simulationFunction(dataFrames[[i]][j,])
            traj[[thisRecord]]$fileNum <- i
        }
    }

    return(traj)
}
