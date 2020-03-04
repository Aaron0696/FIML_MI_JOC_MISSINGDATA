#==== EMB-WLSMV specification ====#
# Model to be fed into lavaan-related commands
lavaanModel <-
  ' F1 =~ ord1*o1 + ord2*o2 + ord3*o3 + cts1*c1 + cts2*c2 + cts3*c3'

# Function to analyze the simulated data
analyze.Data <- function(simDataMissing, numCats)
{
  # Amelia appears to not impute rows with all missing responses,
  # leaving it blank and causing cfa.mi to produce an error.
  # Shall remove the rows with all missing response.
  simDataMissing <- subset(
    simDataMissing,
    is.na(simDataMissing$o1) == FALSE |
      is.na(simDataMissing$o2) == FALSE |
      is.na(simDataMissing$o3) == FALSE |
      is.na(simDataMissing$c1) == FALSE |
      is.na(simDataMissing$c2) == FALSE |
      is.na(simDataMissing$c3) == FALSE
  )
  
  ameliaImputation <- amelia(simDataMissing,
                             m = imputationMI,
                             ords = c("o1", "o2", "o3"))
  # Collating imputed datasets
  ameliaImp <- ameliaImputation$imputations
  
  # Run lavaan with imputed data using runMI
  outputAmelia <- cfa.mi(
    model = lavaanModel,
    data = ameliaImp,
    std.lv = TRUE,
    estimator = "WLSMV",
    parameterization = "theta"
  )
  return(outputAmelia)
}
