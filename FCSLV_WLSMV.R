# Model to be fed into lavaan-related commands
lavaanModel <-
  ' F1 =~ ord1*o1 + ord2*o2 + ord3*o3 + cts1*c1 + cts2*c2 + cts3*c3'

analyze.Data <- function(completeData, numCat)
{
  # container for the results of each imputation
  containimpute <- data.frame(matrix(999, nrow = imputationMI, ncol = 13))
  myimp <- NULL
  # Collating imputed datasets
  for (i in 1:imputationMI)
  {
    myimp[[i]] <- subset(completeData, completeData$IMP == i)
  }
  # run lavaan with imputed data using runMI
  output <- cfa.mi(
    model = lavaanModel,
    data = myimp,
    std.lv = TRUE,
    estimator = "WLSMV",
    parameterization = "theta"
  )
  return(output)
}