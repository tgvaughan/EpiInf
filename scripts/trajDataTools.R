## Scripts useful for plotting and manipulating trajectory data
## using ggplot and friends.

library(dplyr)
library(readr)
library(stringr)

parseTrajectory <- function(trajStr) {
  strValues <- str_split(str_split(trajStr, ",")[[1]], ":", simplify = TRUE)
  res <- list(age = as.numeric(strValues[,1]),
              S = as.numeric(strValues[,2]),
              I = as.numeric(strValues[,3]),
              R = as.numeric(strValues[,4]),
              leap = strValues[,5] == "TL",
              incidence = as.numeric(strValues[,6]),
              Re = as.numeric(strValues[,7]))

  if (dim(strValues)[2]>7)
      res$cumulativeInfections <- as.numeric(strValues[,8])

  return(res)
}

loadTrajectories <- function(filenames, burninFrac=0.1) {
    df <- NULL

    for (filename in filenames) {
        cat(paste("Loading", filename,"..."))
        df_in <- read_tsv(filename, col_types="ic")

        N <- dim(df_in)[1]
        df_in <- df_in[-(1:ceiling(burninFrac*N)),]
        
        for (row in 1:(dim(df_in)[1])) {
            trajStr <- df_in[row,2]
            trajStates <- parseTrajectory(trajStr)
            df <- bind_rows(df,
                            tibble(traj=row,
                                   age=trajStates$age,
                                   S=trajStates$S,
                                   I=trajStates$I,
                                   R=trajStates$R,
                                   leap=trajStates$leap,
                                   incidence=trajStates$incidence,
                                   Re=trajStates$Re,
                                   cumulativeInfections=trajStates$cumulativeInfections))
        }
    }
    
    cat("done.\n")
    
    return(df)
}

gridTrajectories <- function(trajdf, ages) {
    df_grid <- NULL

    my_max <- function(x) {
        if (length(x)==0)
            return(0)
        else
            return (max(x))
    }

    for (grid_age in ages) {

        age_summary <- trajdf %>%
            group_by(traj) %>%
            summarize(
                I=my_max(I[age>grid_age]),
                cumulativeInfections=my_max(cumulativeInfections[age>grid_age]),
                .groups = "drop_last")

        age_summary$age <- grid_age
        df_grid <- bind_rows(df_grid, age_summary)
    }

    return(df_grid)
}
