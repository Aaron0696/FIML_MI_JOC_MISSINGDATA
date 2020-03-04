#==== EMB-MLM specification ====#
# Model to be fed into lavaan-related commands
lavaanModel <-
  ' F1 =~ ord1*o1 + ord2*o2 + ord3*o3 + cts1*c1 + cts2*c2 + cts3*c3'

# Function to analyze the simulated data
analyze.Data <- function(simDataMissing, numCats)
{
  # Coercing the variables back into continuous for Amelia
  simDataMissing$o1 <- as.numeric(simDataMissing$o1)
  simDataMissing$o2 <- as.numeric(simDataMissing$o2)
  simDataMissing$o3 <- as.numeric(simDataMissing$o3)
  
  ameliaImputation <- amelia(simDataMissing,
                             m = imputationMI)
  # Collating imputed datasets
  ameliaImp <- ameliaImputation$imputations
  
  # Run lavaan with imputed data using runMI
  outputAmelia <- cfa.mi(
    model = lavaanModel,
    data = ameliaImp,
    std.lv = TRUE,
    std.ov = TRUE,
    estimator = "MLM",
    parameterization = "theta"
  )
  #fit <- list(lavaan = outputAmelia)
  return(outputAmelia)
}