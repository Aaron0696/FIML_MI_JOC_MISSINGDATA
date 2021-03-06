#==== FCS-WLSMV specification ====#
# Model to be fed into lavaan-related commands
lavaanModel <-
  ' F1 =~ ord1*o1 + ord2*o2 + ord3*o3 + cts1*c1 + cts2*c2 + cts3*c3'

# Function to analyze the simulated data
analyze.Data <- function(simDataMissing, numCat)
{
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
    estimator = "WLSMV",
    parameterization = "theta"
  )
  return(outputMICE)
}