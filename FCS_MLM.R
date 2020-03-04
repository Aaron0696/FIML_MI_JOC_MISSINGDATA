#==== FCS-MLM specification ====#
# Model to be fed into lavaan-related commands
lavaanModel <-
  ' F1 =~ ord1*o1 + ord2*o2 + ord3*o3 + cts1*c1 + cts2*c2 + cts3*c3'

# Function to analyze the simulated data
analyze.Data <- function(simDataMissing, numCat)
{
  # Coercing the variables back into continuous
  simDataMissing$o1 <- as.numeric(simDataMissing$o1)
  simDataMissing$o2 <- as.numeric(simDataMissing$o2)
  simDataMissing$o3 <- as.numeric(simDataMissing$o3)
  
  # Analyzing the simulated dataset with MICE
  # Multiple imputation with the default
  miceImputation <- mice(simDataMissing,
                         m = imputationMI)
  miceImp <- NULL
  
  # Collating imputed datasets
  for (i in 1:imputationMI)
    miceImp[[i]] <- complete(miceImputation,
                             action = i,
                             inc = FALSE)
  
  # Run lavaan with imputed data using runMI
  outputMICE <- cfa.mi(
    model = lavaanModel,
    data = miceImp,
    std.lv = TRUE,
    std.ov = TRUE,
    estimator = "MLM",
    parameterization = "theta"
  )
  return(outputMICE)
}