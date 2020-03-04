# Function to extract relevant data from output files
extract <- function(segment)
{
  df <- data.frame()
  # loop through each core
  for (core in 1:numCore)
  {
    if (length(segment[[core]]) == 0)
    {
      print(paste("Error in Condition ",
                  condition,
                  ", Core ",
                  core,
                  sep = ""))
    } else
    {
      # loop through the lavaan objects within each core
      for (a in seq(1, length(segment[[core]]), by = 2))
      {
        # extract the factor loadings and the SEs
        # create vectors to hold the relevant numbers and names
        estdf <- segment[[core]][[a]][["parameters"]]
        estVec <-
          as.numeric(estdf[estdf$matrix != "Thresholds", ]$Estimate)
        if (length(estVec) != 12)
        {
          print("Error!")
        }
        # standard errors
        seVec <-
          as.numeric(estdf[estdf$matrix != "Thresholds", ]$Std.Error)
        # for column names
        nameVec <- estdf[estdf$matrix != "Thresholds", ]$name
        # whether the SEs are suspect
        sesusVec <- segment[[core]][[a]][["seSuspect"]][1:12]
        # whether the covariance matrix is info definite
        infDef <- segment[[core]][[a]][["infoDefinite"]]
        # degree of freedom
        ChiDofVec <- segment[[core]][[a]][["ChiDoF"]]
        # chi square value
        ChiVec <- segment[[core]][[a]][["Chi"]]
        # p value
        pVec <- segment[[core]][[a]][["p"]]
        # CFI, TFI, RMSEA
        cfiVec <- segment[[core]][[a]][["CFI"]]
        tliVec <- segment[[core]][[a]][["TLI"]]
        rmseaLVec <-
          as.numeric(segment[[core]][[a]][["RMSEACI"]][["lower"]])
        rmseaHVec <-
          as.numeric(segment[[core]][[a]][["RMSEACI"]][["upper"]])
        rmseaVec <- as.numeric(segment[[core]][[a]][["RMSEA"]])
        # status of the estimator, indicates if convergence was reached
        statusVec <-
          as.character(segment[[core]][[a]][["statusCode"]])
        
        # making column names for the data.frame
        senameVec <- paste(nameVec,
                           "SE",
                           sep = ".")
        sesusnameVec <- paste(nameVec,
                              "SEsus",
                              sep = ".")
        
        # bind all into one row
        allstats <- c(
          estVec,
          seVec,
          sesusVec,
          infDef,
          ChiDofVec,
          ChiVec,
          pVec,
          cfiVec,
          tliVec,
          rmseaLVec,
          rmseaHVec,
          rmseaVec,
          statusVec
        )
        
        # convert to dataframe and name the columns
        allstats <-
          data.frame(t(allstats), stringsAsFactors = FALSE)
        names(allstats) <- c(
          nameVec,
          senameVec,
          sesusnameVec,
          "infDef",
          "ChiDof",
          "Chi",
          "pChi",
          "cfi",
          "tli",
          "rmseaL",
          "rmseaH",
          "rmsea",
          "status"
        )
        # bind it to the empty container "df"
        df <-
          rbind(df,
                allstats,
                make.row.names = FALSE,
                row.names = NULL)
        # making sure column names are there
        names(df) <- c(
          nameVec,
          senameVec,
          sesusnameVec,
          "infDef",
          "ChiDof",
          "Chi",
          "pChi",
          "cfi",
          "tli",
          "rmseaL",
          "rmseaH",
          "rmsea",
          "status"
        )
        
        paste("Core ",
              core,
              " Number ",
              a,
              " Done",
              sep = "")
      }
    }
  }
  return(df)
}
