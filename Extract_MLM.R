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
      for (a in seq(1,length(segment[[core]]), by = 2))
      {
        # save a summary of the lavaan results
        summ <- summary(segment[[core]][[a]], standardized = TRUE)
        # factor loadings
        loadVec <- summ$est[grep("=~", summ$op)]
        # error variances
        errorVec <-  summ$est[grep("~~", summ$op)]
        # exclude the error variance of the factor
        errorVec <- errorVec[-length(errorVec)]
        # standard error
        seVec <- summ$se[grep("~~|=~", summ$op)]
        # exclude se of the error variance of the factor
        seVec <- seVec[-length(seVec)]
        
        # name vectors
        nameVec <- segment[[core]][[a]]@Model@dimNames[[2]][[1]]
        errornameVec <- paste(nameVec,
                              "error", 
                              sep = ".")
        senameVec <- paste(nameVec,
                           "SE", 
                           sep = ".")
        errorsenameVec <- paste(nameVec,
                                "errorSE",
                                sep = ".")
        
        temp <- unlist(segment[[core]][[a]]@convergence)
        # Number of Convergence
        tempConverge <- temp[seq(from = 1, to = imputationMI*4, by = 4)]
        numConverge <- sum(tempConverge, na.rm = TRUE)
        # Number of SE = TRUE
        tempSE <- temp[seq(from = 2, to = imputationMI*4, by = 4)]
        numSE <- sum(tempSE, na.rm = TRUE)
        # Number of Heywood.lv
        tempHeylv <- temp[seq(from = 3, to = imputationMI*4, by = 4)]
        numHeylv <- sum(tempHeylv, na.rm = TRUE)
        # Number of Heywood.ov
        tempHeyov <- temp[seq(from = 4, to = imputationMI*4, by = 4)]
        numHeyov <- sum(tempHeyov, na.rm = TRUE)
        # All Convergence/Heywood counts in one row
        oneRow <- c(numConverge, numSE, numHeylv, numHeyov)
        # Name for the above vectors
        nameConvVec  <- c("Convergence", "SETrue", "Heywood.lv", "Heywood.ov")
        
        mengrubin <- tryCatch({
          lavTestLRT.mi(segment[[core]][[a]],
                        asymptotic = TRUE, 
                        test = "D3")
        }, error = function(e){
          rep(NA,16)
        })
        
        # in the case that the model fit is perfect
        if (length(mengrubin) == 12)
        {
          mengrubin <- rep("Perfect",16)
        }
        
        
        names(mengrubin) <- paste(c("chisq",
                                    "df",
                                    "p",
                                    "ariv",
                                    "fmi",
                                    "npar",
                                    "ntotal",
                                    "logl",
                                    "unrestricted.logl",
                                    "aic",
                                    "bic",
                                    "bic2",
                                    "chisq.scaled",
                                    "df.scaled",
                                    "p.scaled",
                                    "chisq.scale.factor"),
                                  "mengrubin", sep = ".")
        
        Lirobust <- tryCatch({
          lavTestLRT.mi(segment[[core]][[a]],
                        asymptotic = TRUE, 
                        test = "D2",
                        pool.robust = TRUE)
        }, error = function(e){
          rep(NA,12)
        })
        
        names(Lirobust) <- paste(c("chisq",
                                   "df",
                                   "p",
                                   "ariv",
                                   "fmi",
                                   "npar",
                                   "ntotal",
                                   "chisq.scaled",
                                   "df.scaled",
                                   "p.scaled",
                                   "ariv.scaled",
                                   "fmi.scaled"),
                                 "lirobust", sep = ".")
        print(core)
        print(a)
        Li <- tryCatch({
          lavTestLRT.mi(segment[[core]][[a]],
                        asymptotic = TRUE, 
                        test = "D2", 
                        pool.robust = FALSE)
        }, error = function(e){
          rep(NA,11)
        })
        
        # in the case that the model fit is perfect
        if (length(Li) == 7) 
        {
          Li <- rep("Perfect", 11)
        }
        
        names(Li) <- paste(c("chisq",
                             "df",
                             "p",
                             "ariv",
                             "fmi",
                             "npar",
                             "ntotal",
                             "chisq.scaled",
                             "df.scaled",
                             "p.scaled",
                             "chisq.scaling.factor"),
                           "li", sep = ".")
        
        # Bind them all vectors into one row
        allstats <- c(loadVec,
                      errorVec,
                      seVec,
                      oneRow,
                      mengrubin,
                      Li,
                      Lirobust)
        
        allstats <- data.frame(t(allstats), 
                               stringsAsFactors = FALSE)
        names(allstats) <- c(nameVec,
                             errornameVec,
                             senameVec,
                             errorsenameVec,
                             nameConvVec,
                             names(mengrubin),
                             names(Li),
                             names(Lirobust))
        df <- rbind(df, 
                    allstats, 
                    make.row.names = FALSE, 
                    row.names = NULL)
        names(df) <- c(nameVec,
                       errornameVec,
                       senameVec,
                       errorsenameVec,
                       nameConvVec,
                       names(mengrubin),
                       names(Li),
                       names(Lirobust))
        
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
