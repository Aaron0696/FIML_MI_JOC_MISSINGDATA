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

# Function to analyze the simulated data
analyze.Data <- function(simDataMissing, numCat)
{
  # OpenMX requires the setting of thresholds for analysis with ordinal variables
  # nThresh is the number of levels in the ordinal variable minus 1
  thresholds <- mxThreshold(
    vars = c("o1", "o2", "o3"),
    nThresh = c(numCat - 1),
    free = c(TRUE)
  )
  
  # Analyzing the simulated dataset with OpenMx
  dataRaw <- mxData(observed = simDataMissing,
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
  factorFit <- mxRun(OMXmodel)
  
  # Check the status code if it is either 0 or 1, run mxTryHard if it is not.
  if (!(factorFit$output$status[[1]] %in% c(0, 1))) {
    factorFit <- mxTryHardOrdinal(factorFit, extraTries = 30)
  }
  # Check again. It still does not work, output NA
  if (!(factorFit$output$status[[1]] %in% c(0, 1))) {
    cat("Output NA As Results")
  } else {
    cat("Get The Results")
  }
  satFit <- omxSaturatedModel(OMXmodel, run = TRUE)
  fit <- summary(factorFit, refModels = satFit)
  return(fit)
}