#==== JOC-FIML model specification ====#
# Residual Variances
resVars <- mxPath(
  from = c("o1", "o2", "o3", "c1", "c2", "c3"),
  arrows = 2,
  free = c(FALSE, FALSE, FALSE, TRUE, TRUE, TRUE),
  values = c(1, 1, 1, 1, 1, 1),
  labels = c("Eo1", "Eo2", "Eo3", "Ec1", "Ec2", "Ec3")
)

# Means
means <- mxPath(
  from = "one",
  to = c("o1", "o2", "o3", "c1", "c2", "c3", "F1"),
  arrows = 1,
  free = c(FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, FALSE),
  values = c(0),
  labels = c(
    "mean_o1",
    "mean_o2",
    "mean_o3",
    "mean_c1",
    "mean_c2",
    "mean_c3",
    "mean_F1"
  )
)

# Latent Variances
latVars <- mxPath(
  from = c("F1"),
  arrows = 2,
  connect = "unique.pairs",
  free = FALSE,
  values = c(1),
  labels = c("varF1")
)

# Freed all loadings
# Factor loadings for all variables
facLoadsX16 <- mxPath(
  from = "F1",
  to = c("o1", "o2", "o3", "c1", "c2", "c3"),
  arrows = 1,
  free = c(TRUE, TRUE, TRUE, TRUE, TRUE, TRUE),
  values = c(0.3, 0.5, 0.7, 0.7, 0.8, 0.85),
  labels = c("Lo1", "Lo2", "Lo3", "Lc1", "Lc2", "Lc3")
)

analyze.Data <- function(completeData, numCat)
{
  # OpenMX requires the setting of thresholds for analysis with ordinal variables
  # nThresh is the number of levels in the ordinal variable minus 1
  # trueThreshold[c(-1,-numCat-1)]
  thresholds <- mxThreshold(
    vars = c("o1", "o2", "o3"),
    nThresh = c(numCat - 1),
    free = c(TRUE)
  )
  
  # container for the results of each imputation
  containimpute <- data.frame(matrix(999, nrow = imputationMI, ncol = 13))
  # collating imputed datasets, for each imputed dataset...
  for (i in 1:imputationMI)
  {
    # extract complete data
    mydata <- completeData[completeData$IMP == i,]
    mydata[,c("o1","o2","o3")] <- lapply(mydata[,c("o1","o2","o3")], factor, level = 1:numCat, ordered = TRUE)
    
    # analyzing the simulated dataset with OpenMx
    dataRaw <- mxData(observed = mydata,
                      type = "raw")
    OMXmodel <- mxModel(
      "Simulated",
      type = "RAM",
      manifestVars = c("o1", "o2", "o3", "c1", "c2", "c3"),
      latentVars = c("F1"),
      dataRaw,
      resVars,
      latVars,
      facLoadsX16,
      means,
      thresholds
    )
    # fit the model
    factorFit <- mxRun(OMXmodel)
    ## Check the status code if it is either 0 or 1, run mxTryHard if it is not.
    if (!(factorFit$output$status[[1]] %in% c(0, 1))) {
      factorFit <- mxTryHardOrdinal(factorFit, extraTries = 30)
    }
    ## Check again. It still does not work, output NA
    if (!(factorFit$output$status[[1]] %in% c(0, 1))) {
      cat("Output NA As Results")
    } else {
      cat("Get The Results")
    }
    
    # get factor loading and SE from summary()
    fit <- summary(factorFit)
    # extract the parameters as a dataframe saved in params
    params <- fit[["parameters"]]
    # take the factor loadings estimates (contains "L" in the name) and insert them into the first six columns
    containimpute[i,1:6] <- t(params[grep("L" ,params$name), c("Estimate")])
    # do the same for SE, put them in the next six columns (12 in total)
    containimpute[i,7:12] <- t(params[grep("L" ,params$name), c("Std.Error")])
    # use the last column to note if convergence is reached
    # convergence is reached 
    containimpute[i,13] <- factorFit$output$status[[1]] %in% c(0, 1)
  }
  
  # rename the columns
  names(containimpute) <- c("o1", "o2", "o3", "c1", "c2", "c3","o1.SE", "o2.SE", "o3.SE", "c1.SE", "c2.SE", "c3.SE", "Convergence")
  
  # if there is NA in containimpute, change convergence to FALSE
  containimpute[["Convergence"]] <- ifelse(rowSums(is.na(containimpute)) > 0, 
                                           FALSE, 
                                           containimpute[["Convergence"]])
  # remove the imputations which did not converge
  containimpute <- containimpute[containimpute$Convergence != 0,]
  
  # output is the vector containing the pooled estimate and SE for each variable
  output <- c(rep(999,13))
  names(output) <- paste0(names(containimpute), ".pooled")
  # create a vector of variables to loop through
  myvec <- c("o1","o2","o3","c1","c2","c3")
  
  for(r in 1:length(myvec))
  {
    pooled <- pool.scalar(Q = containimpute[[myvec[r]]], # Q: A vector of univariate estimates of m repeated complete data analyses.
                          U = containimpute[[paste0(myvec[r], ".SE")]]) # U: A vector containing the corresponding m variances of the univariate estimates.
    output[r] <- pooled$qbar
    output[r + 6] <- pooled$t
  }
  
  # count number of convergence
  output[13] <- sum(containimpute$Convergence) 
  return(output)
}