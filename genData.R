# function to generate data with the requested traits
gen.Data <-
  function(sampleSize,
           missingMech,
           missingProp,
           numCat,
           aSym)
  {
    # simulating the raw and continuous data
    simDataMissing <- simulateData(
      simulationModel,
      sample.nobs = sampleSize,
      return.type = "data.frame",
      return.fit = TRUE,
      model.type = "cfa"
    )
    
    # getting true thresholds for each combination
    if (numCat == 2 & aSym == 1)
    {
      trueThreshold <- symThreshold2
    } else if (numCat == 2 & aSym == 2)
    {
      trueThreshold <- MsymThreshold2
    } else if (numCat == 2 & aSym == 3)
    {
      trueThreshold <- SsymThreshold2
      
      ####### 3 Categories #######
    } else if (numCat == 3 & aSym == 1)
    {
      trueThreshold <- symThreshold3
    } else if (numCat == 3 & aSym == 2)
    {
      trueThreshold <- MsymThreshold3
    } else if (numCat == 3 & aSym == 3)
    {
      trueThreshold <- SsymThreshold3
      
      ####### 5 Categories #######
    } else if (numCat == 5 & aSym == 1)
    {
      trueThreshold <- symThreshold5
    } else if (numCat == 5 & aSym == 2)
    {
      trueThreshold <- MsymThreshold5
    } else if (numCat == 5 & aSym == 3)
    {
      trueThreshold <- SsymThreshold5
      
      ####### 7 Categories #######
    } else if (numCat == 7 & aSym == 1)
    {
      trueThreshold <- symThreshold7
    } else if (numCat == 7 & aSym == 2)
    {
      trueThreshold <- MsymThreshold7
    } else if (numCat == 7 & aSym == 3)
    {
      trueThreshold <- SsymThreshold7
    } else
    {
      print("ERROR IN THRESHOLD ASSIGNMENT")
    }
    
    # create missing data
    # create random numbers to use for MCAR missingness
    o1rand <- runif(sampleSize)
    o2rand <- runif(sampleSize)
    o3rand <- runif(sampleSize)
    c1rand <- runif(sampleSize)
    c2rand <- runif(sampleSize)
    c3rand <- runif(sampleSize)
    
    # case is missing if its random value is less than missing proportion
    simDataMissing$o1[o1rand < missingProp] <- NA
    simDataMissing$o2[o2rand < missingProp] <- NA
    simDataMissing$o3[o3rand < missingProp] <- NA
    simDataMissing$c1[c1rand < missingProp] <- NA
    simDataMissing$c2[c2rand < missingProp] <- NA
    simDataMissing$c3[c3rand < missingProp] <- NA
    
    # cut into ordinal variables
    simDataMissing[, 1] <-
      cut(
        simDataMissing[, 1],
        breaks = c(-Inf, trueThreshold, Inf),
        labels = c(1:numCat),
        right = TRUE
      )
    simDataMissing[, 2] <-
      cut(
        simDataMissing[, 2],
        breaks = c(-Inf, trueThreshold, Inf),
        labels = c(1:numCat),
        right = TRUE
      )
    simDataMissing[, 3] <-
      cut(
        simDataMissing[, 3],
        breaks = c(-Inf, trueThreshold, Inf),
        labels = c(1:numCat),
        right = TRUE
      )
    
    # coercing the variables into ordered variables
    simDataMissing[, 1] <- as.ordered(simDataMissing[, 1])
    simDataMissing[, 2] <- as.ordered(simDataMissing[, 2])
    simDataMissing[, 3] <- as.ordered(simDataMissing[, 3])
    return(simDataMissing)
  }