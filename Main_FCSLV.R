# packages
library(semTools)
library(mice)
library(Amelia)
library(OpenMx)
library(foreach)
library(doParallel)
library(doRNG)

# Setting Approach, Number Of Replications and Cores -----------------------------------------

# approach to use
approach <- "FCSLV_WLSMV"
# approach <- "FCSLV_FIML"

# number of cores to use
numCore <- 3
# number of replications in each condition, in each core, e.g if core = 4 and itert = 3, 12 replications will be ran for each condition
itert <- 4
numCore <- itert # trick the extract function to loop through each iteration rather than each core for FCSLV-WLSMV
registerDoParallel(cores = numCore)

# for replicability
set.seed(123456)
# load functions used for data generation
source("genData.R")

# obtained from "Multivariate normal maximum
# likelihood with both ordinal and continuous
# variables, and data missing at random". Pritikin, Brick & Neale (2018)
mxOption(NULL, "Default optimizer", "SLSQP")  # CSOLNP has trouble with thresholds

# Set Manipulated Dataset Factors -----------------------------------------
# sample size
sampleSize <- c(500)
# mechanism of missingness ["MCAR", "MAR"]
missingMech <- c("MCAR")
# proportion of missing data
missingProp <- c(0.1)
# number of categories of ordinal data
numCat <- c(2, 3)
# threshold for ordinal variables
# symmetrical
symThreshold2 <- c(0)
symThreshold3 <- c(-0.83, 0.83)
symThreshold5 <- c(-1.50, -0.50, 0.50, 1.50)
symThreshold7 <- c(-1.79, -1.07, -0.36, 0.36, 1.07, 1.79)
# moderate asymmetry
MsymThreshold2 <- c(0.36)
MsymThreshold3 <- c(-0.50, 0.76)
MsymThreshold5 <- c(-0.70, 0.39, 1.16, 2.05)
MsymThreshold7 <- c(-1.43, -0.43, 0.38, 0.94, 1.44, 2.54)
# severe asymmetry
SsymThreshold2 <- c(1.04)
SsymThreshold3 <- c(0.58, 1.13)
SsymThreshold5 <- c(0.05, 0.44, 0.84, 1.34)
SsymThreshold7 <- c(-0.25, 0.13, 0.47, 0.81, 1.18, 1.64)
# number of imputations for Multiple Imputation
imputationMI <- 10
# a dummy to be used to input the degree of asymmetry
# 1 for symmetric
# 2 for moderate
# 3 for severe
aSymDummy <- seq(from = 1, to = 3, by = 1)
# creating the matrix showing all possible conditions
conditionsMatrix <- expand.grid(
  sampleSize = sampleSize,
  missingMech = missingMech,
  missingProp = missingProp,
  aSym = aSymDummy,
  numCat = numCat
)
# Creating columns for the threshold
conditionsMatrix$threshold1 <- NA
conditionsMatrix$threshold2 <- NA
conditionsMatrix$threshold3 <- NA
conditionsMatrix$threshold4 <- NA
conditionsMatrix$threshold5 <- NA
conditionsMatrix$threshold6 <- NA

# fill in threshold values based on numCat and aSym
for (r in 1:nrow(conditionsMatrix))
{
  # 2 Categories
  
  if (conditionsMatrix$numCat[r] == 2 & conditionsMatrix$aSym[r] == 1)
  {
    conditionsMatrix[r, 6] <- symThreshold2
  } else if (conditionsMatrix$numCat[r] == 2 &
             conditionsMatrix$aSym[r] == 2)
  {
    conditionsMatrix[r, 6] <- MsymThreshold2
  } else if (conditionsMatrix$numCat[r] == 2 &
             conditionsMatrix$aSym[r] == 3)
  {
    conditionsMatrix[r, 6] <- SsymThreshold2
    
    # 3 Categories
    
  } else if (conditionsMatrix$numCat[r] == 3 &
             conditionsMatrix$aSym[r] == 1)
  {
    conditionsMatrix[r, 6:7] <- symThreshold3
  } else if (conditionsMatrix$numCat[r] == 3 &
             conditionsMatrix$aSym[r] == 2)
  {
    conditionsMatrix[r, 6:7] <- MsymThreshold3
  } else if (conditionsMatrix$numCat[r] == 3 &
             conditionsMatrix$aSym[r] == 3)
  {
    conditionsMatrix[r, 6:7] <- SsymThreshold3
    
    # 5 Categories
    
  } else if (conditionsMatrix$numCat[r] == 5 &
             conditionsMatrix$aSym[r] == 1)
  {
    conditionsMatrix[r, 6:9] <- symThreshold5
  } else if (conditionsMatrix$numCat[r] == 5 &
             conditionsMatrix$aSym[r] == 2)
  {
    conditionsMatrix[r, 6:9] <- MsymThreshold5
  } else if (conditionsMatrix$numCat[r] == 5 &
             conditionsMatrix$aSym[r] == 3)
  {
    conditionsMatrix[r, 6:9] <- SsymThreshold5
    
    # 7 Categories
    
  } else if (conditionsMatrix$numCat[r] == 7 &
             conditionsMatrix$aSym[r] == 1)
  {
    conditionsMatrix[r, 6:11] <- symThreshold7
  } else if (conditionsMatrix$numCat[r] == 7 &
             conditionsMatrix$aSym[r] == 2)
  {
    conditionsMatrix[r, 6:11] <- MsymThreshold7
  } else if (conditionsMatrix$numCat[r] == 7 &
             conditionsMatrix$aSym[r] == 3)
  {
    conditionsMatrix[r, 6:11] <- SsymThreshold7
  } else
  {
    print("ERROR IN THRESHOLD ASSIGNMENT")
  }
}

# population model
simulationModel <-
  'F1 =~ 0.3*o1 + 0.5*o2 + 0.7*o3 + 0.7*c1 + 0.8*c2 + 0.85*c3
  # Fix variance of latent factor to 1
  F1 ~~ 1*F1

  # Fixing error variances
  o1 ~~ 1*o1
  o2 ~~ 1*o2
  o3 ~~ 1*o3
  c1 ~~ 0.51*c1
  c2 ~~ 0.36*c2
  c3 ~~ 0.2775*c3'


# Generate Data -----------------------------------------------------------------
# create .dat files to be used with blimp software
# .dat files should appear in the working directory
for(r in 1:nrow(conditionsMatrix))
{
  counter <- 1
  replicate(itert, expr = {
    mydata <- gen.Data(sampleSize = conditionsMatrix$sampleSize[r],
                       missingMech = "MCAR",
                       missingProp = conditionsMatrix$missingProp[r],
                       numCat = conditionsMatrix$numCat[r],
                       aSym = conditionsMatrix$aSym[r])
    mydata <- data.frame(lapply(mydata, as.numeric))
    mydata[is.na(mydata)] <- 999
    write.table(mydata, file = paste0(r,"_","data",counter,".dat"), sep = " ", col.names = FALSE, row.names = FALSE)
    counter <<- counter + 1
  })
}


# Generate blimp instructions ---------------------------------------------
# automatically generate the instructions for blimp in .imp files
contain <- c() # container to store text for .bat file
for(r in 1:nrow(conditionsMatrix))
{
  for(i in 1:itert)
  {
    fileConn<-file(paste0(r,"_imp",i,".imp"))
    writeLines(c(paste0("DATA: ",r,"_data",i,".dat;"),
                 "VARIABLES: o1 o2 o3 c1 c2 c3;",
                 "NOMINAL: ;",
                 "ORDINAL: o1 o2 o3;",
                 "MISSING: 999;",
                 "FCS: o1 o2 o3 c1 c2 c3;",
                 "SEED: 123456;",
                 "NIMPS: 50;",
                 "BURN: 5000;",
                 "THIN: 2500;",
                 "CHAINS: 2 processors 2;",
                 paste0("OUTFILE: ", r, "_myimps",i,".csv;"),
                 "OPTIONS: psr latent;"
    ), 
    fileConn)
    close(fileConn)
    contain <- c(contain,paste0("blimp ",getwd(),"/",r,"_imp",i,".imp"))
  }
}

# Generate .bat file to run all .imp files --------------------------------
# write .bat file
fileConn <- file("Blimp.bat")
writeLines(c(contain), 
           fileConn)
close(fileConn)

# Run blimp to generate imputed data -------------------------------------
#! done outside of R
# Load imputed data and analyze ------------------------------------------
source(paste0(approach, ".R")) # load analyze.Data() function

for (condition in 1:nrow(conditionsMatrix))
{
  name <- paste(approach, "_output", condition, sep = "")
  mylist <- list() # list to contain the outputs of each iteration
  for (r in 1:itert)
  {
    # load data in
    data <- read.csv(paste0(condition,"_myimps",r,".csv"))
    names(data) <- c("IMP","o1","o2","o3","c1","c2","c3")
    fitInfo <- list(fit = tryCatch({
      analyze.Data(data, conditionsMatrix$numCat[condition])},
      error = function(e) {
        "Error"
      }))
    mylist[[r]] <- fitInfo
    # save output to an RDS file
    saveRDS(mylist, paste(name,
                          "RDS",
                          sep = "."))
  }
}

# Extract data if FCSLV_WLSMV is used -------------------------------------
if(approach == "FCSLV_WLSMV")
{
  source("Extract_WLSMV.R")
  # summarize the results from each output file into a
  # dataframe with each row representing one
  # replication in the specific condition
  for(num in 1:nrow(conditionsMatrix)){
    # create the name of the output
    outname <- paste(approach, "_output", num, sep = "")
    # read the file in and assign it to the segment object
    segment <- readRDS(paste(outname, ".RDS", sep = ""))
    df <- extract(segment)
    # add in columns detailing the factors used to generate
    # the data
    df$sampleSize <- conditionsMatrix$sampleSize[num]
    df$missingProp <- conditionsMatrix$missingProp[num]
    df$numCat <- conditionsMatrix$numCat[num]
    df$aSym <- conditionsMatrix$aSym[num]
    df$approach <- approach
    # save df into the desktop
    saveRDS(df, 
            paste(approach, 
                  "_Summarize_",
                  num,
                  ".RDS",
                  sep = ""))
  }
}


