# Libraries ---------------------------------------------------------------
# for latent variable analysis
library(semTools)
# for multiple imputation using FCS
library(mice)
# for multiple imputation using EMB
library(Amelia)
# for JOC-FIML
library(OpenMx)
# for data wrangling
library(reshape2)
library(dplyr)
# for visualization
library(ggplot2)
library(corrplot)
# for loading XPT files
library(Hmisc)

# Loading Data ------------------------------------------------------------
# Physical Functioning Questionnaire
# 2013 to 2014 cycle
PFQdat <- sasxport.get("https://wwwn.cdc.gov/Nchs/Nhanes/2013-2014/PFQ_H.XPT")

# Grip Test
# 2013 to 2014 cycle
GRIPdat <- sasxport.get("https://wwwn.cdc.gov/Nchs/Nhanes/2013-2014/MGX_H.XPT")

# Demographic Data
# 2013 to 2014 cycle
DEMOdat <- sasxport.get("https://wwwn.cdc.gov/Nchs/Nhanes/2013-2014/DEMO_H.XPT")

# merge dataframes together using full outer-join to 
# obtain `mydata`, a consolidated dataset for analysis.

# full outer join these two datasets by ID number
mydata <- merge(GRIPdat,
                PFQdat,
                by = "seqn", 
                all = "true")
# perform the same outer join with the joint data from above and the demographic data
mydata <- merge(mydata,
                DEMOdat,
                by = "seqn", 
                all = "true")

# only respondents who responded 1/YES to either pfq049, 
# pfq051, pfq054, pfq057 were shown the pfq061 questions

# respondents who had missing data on all pfq061 questions 
# responded 2/NO to pfq049, pfq051, pfq054 and pfq057
# checkdat was created to check if this is true by
# only including those with missing responses on 
# pfq061a
checkdat <- subset(PFQdat,
                   is.na(PFQdat$pfq061a))
# we can see that the rest of the pfq061 questions were completely 
# missing and there were only "2" responses to  
# pfq049, pfq051, pfq054 and pfq057
summary(checkdat)


# Creating Dataframe for Analysis -----------------------------------------
# extract only variables that we are interested in (denoted by *)
# mgxh1t1 - Grip strength of hand1 at test1
# mgxh1t1e - Whether the participant exerted maximal(1) or questionable effort(2) for hand1 at test1
# mgxh1t2 - Grip strength of hand1 at test2
# mgxh1t2e - Whether the participant exerted maximal(1) or questionable effort(2) for hand1 at test2
# mgxh1t3 - Grip strength of hand1 at test3
# mgxh1t3e - Whether the participant exerted maximal(1) or questionable effort(2) for hand1 at test3

# mgxh2t1 - Grip strength of hand2 at test1
# mgxh2t1e - Whether the participant exerted maximal(1) or
# mgxh2t2 - Grip strength of hand2 at test1
# mgxh2t2e - Whether the participant exerted maximal(1) or questionable effort(2) for hand2 at test2
# mgxh2t3 - Grip strength of hand2 at test3
# mgxh2t3e - Whether the participant exerted maximal(1) or questionable effort(2) for hand2 at test3

# *mgdcgsz - Combined grip strength that is the sum of the highest grip strength recorded for each hand

# *pfq061e - Lifting or carrying difficulty
# pfq061k - Using fork, knife, drinking from cup
# pfq061l - Dressing yourself difficulty
# pfq061o - Reaching up over head difficulty
# *pfq061p - Grasp/holding small objects difficulty
# *pfq061t - Push or pull large objects difficulty
# ridageyr - Age at time of examination in years

# this command reads:
# retrieve the columns of mydata that contains the characters
# ridageyr or mgdcgsz or pfq061e or pfq061p or pfq061t
mydata <- data.frame(mydata[,grep("ridageyr|mgdcgsz|pfq061[e,p,t]",
                                  names(mydata))])
# convert the questionaire items into factors using a for loop
# for each column in mydata
for(r in 1:ncol(mydata)){
  
  # if the column possess less than 10 unique values
  if(length(unique(mydata[,r])) < 10){
    
    # then convert it to a factor
    mydata[,r] <- factor(mydata[,r],
                         # only take 1, 2, 3 and 4 as valid categories
                         # other categories detected will be automatically classified as NA
                         # the levels argument will need to be updated for other datasets
                         # with different number of categories
                         levels = c(1,2,3,4),
                         ordered = TRUE)
    
  }
}
# renaming columns to make data more user-friendly
names(mydata) <- c("CGrip",
                   "LiftD",
                   "GraspD",
                   "PshPlD",
                   "Age")
# delete cases where participants were less than 20 years old
mydata <- mydata[mydata$Age > 20,]

# delete cases where participants had no responses on Cgrip, LiftD, GraspD and PshPlD
mydata <- subset(mydata,
                 !is.na(mydata$CGrip) | !is.na(mydata$LiftD) | !is.na(mydata$GraspD) | !is.na(mydata$PshPlD))

# add an ID variable that is an index of the row number
mydata$ID <- 1:nrow(mydata)

summary(mydata)

# creating another dataframe with continuous versions
# of the scale variables
mydata.cts <- mydata
mydata.cts$LiftD <- as.numeric(as.character(mydata.cts$LiftD))
mydata.cts$GraspD <- as.numeric(as.character(mydata.cts$GraspD))
mydata.cts$PshPlD <- as.numeric(as.character(mydata.cts$PshPlD))

# Univariate Distribution -------------------------------------------------
melted <- melt(mydata, id.vars = "ID")
melted$value <- as.numeric(melted$value)

# plotting of histograms for the questionnaire items
ggplot(data = melted[grep("PshPl|Grasp|Lift", melted$variable),], aes(x = value, fill = Missing)) + 
  geom_histogram(fill = "steelblue", color = "black", binwidth = 1) + 
  facet_wrap(~variable) + 
  xlab("Response") + 
  ylab("Frequency") + 
  theme(strip.text = element_text(size = 12))
# plotting of histograms for CGrip
ggplot(data = melted[grep("CGrip",melted$variable),], aes(x = value)) + 
  geom_histogram(fill = "steelblue", color = "black", binwidth = 1) + 
  xlab("Combined Grip Strength (KG)") + 
  ylab("Frequency") + 
  facet_wrap(~variable) + 
  theme(strip.text = element_text(size = 12))
# collapse the last two categories into one to form three category responses
# this is due to the small frequencies in the last category
levels(mydata$LiftD) <- c("1","2","3","3")
levels(mydata$GraspD) <- c("1","2","3","3")
levels(mydata$PshPlD) <- c("1","2","3","3")


# Bivariate Distribution --------------------------------------------------
# create a temporary dataframe without ID
forcorr <- mydata[,-grep("^ID", names(mydata))]
# convert all columns to numeric after converting to character
# converting directly from factor to numeric will create
# wrong categories if there are some unobserved values or
# if factors are not numeric or do not begin from 1
forcorr <- apply(forcorr,
                 2,
                 function(x)as.numeric(as.character(x)))

# plot correlation matrix using listwise deletion
corrplot.mixed(cor(na.omit(forcorr)), lower = "number", upper = "ellipse", tl.col = "black", tl.cex = 1, tl.pos = "lt")

# Missing Patterns --------------------------------------------------------
# get the missing data as a matrix
missmatrix <- md.pattern(mydata, rotate.names = TRUE, plot = FALSE)
# rename the last column to the number of missing variables of the specific pattern
names(missmatrix)[ncol(missmatrix)] <- "NumMisVar"
# enter the row number of missmatrix as another column
missmatrix <- cbind(missmatrix,as.numeric(row.names(missmatrix)))
# convert to a dataframe
missmatrix <- as.data.frame(missmatrix, row.names = FALSE)
# rename the last column to count which represents the frequency of the specific pattern
names(missmatrix)[ncol(missmatrix)] <- "Count"
# show top 5 in decreasing order of count
head(missmatrix[order(missmatrix[,grep("Count", names(missmatrix))], decreasing = TRUE),],10)

# FCS-WLSMV ---------------------------------------------------------------
model <- '
PF =~ CGrip + LiftD + GraspD + PshPlD
'
# FCS with mice using only analyzed variables
miceImputation <- mice(mydata[,grep("Lift|Grasp|PshPl|CGrip",
                                    names(mydata))],
                       m = 50,
                       seed = 910313)
# creating an object to store the 50 imputed datasets
miceImp <- NULL

# assessing convergence
plot(miceImputation)
# collating imputed datasets
for (i in 1:50){
  miceImp[[i]] <- complete(miceImputation,
                           action = i,
                           inc = FALSE)
}
# run lavaan with imputed data using runMI
outputMICE <- cfa.mi(
  model = model,
  data = miceImp,
  std.lv = TRUE,
  estimator = "WLSMV",
  parameterization = "theta"
)

# view results
summary(outputMICE)

# get chi-square statistics
lavTestLRT.mi(outputMICE,
              test = "D2",
              asymptotic = TRUE,
              pool.robust = FALSE)

# EMB-WLSMV ---------------------------------------------------------------
model <- '
PF =~ CGrip + LiftD + GraspD + PshPlD
'
# EMB with Amelia using only analyzed variables
ameliaImputation <- amelia(mydata[,grep("Lift|Grasp|PshPl|CGrip",
                                        names(mydata))],
                           m = 50,
                           ords = c("LiftD","GraspD","PshPlD"))
# Collating imputed datasets
ameliaImp <- ameliaImputation$imputations

# run lavaan with imputed data using runMI
outputAmelia <- cfa.mi(
  model = model,
  data = ameliaImp,
  std.lv = TRUE,
  estimator = "WLSMV",
  parameterization = "theta"
)

# view results
summary(outputAmelia)

# get chi-square statistics
lavTestLRT.mi(outputAmelia,
              test = "D2",
              pool.robust = TRUE,
              asymptotic = TRUE)

# FCSLV-WLSMV ---------------------------------------------------------------
model <- '
PF =~ CGrip + LiftD + GraspD + PshPlD
'
# extract data to be fed into blimp
temp <- data.frame(lapply(mydata,as.numeric))
# replace NA with numeric code for missingness
temp[is.na(temp)] <- 999
temp %>% select(-ID) %>%
  write.table(file = "empiricalblimp.csv", sep = ",", col.names = FALSE, row.names = FALSE)

# read imputed dataset
blimp_data <- read.csv("blimpout.csv")
names(blimp_data) <- c("IMP","CGrip","LiftD","GraspD","PshPlD","Age")
blimpimp <- NULL
# collating imputed datasets
for (i in 1:50)
{
  blimpimp[[i]] <- subset(blimp_data, blimp_data$IMP == i)
  # change to ordered
  blimpimp[[i]]$LiftD <- factor(blimpimp[[i]]$LiftD, levels = 1:3, ordered = TRUE)
  blimpimp[[i]]$GraspD <- factor(blimpimp[[i]]$GraspD, levels = 1:3, ordered = TRUE)
  blimpimp[[i]]$PshPlD <- factor(blimpimp[[i]]$PshPlD, levels = 1:3, ordered = TRUE)
}

# run lavaan with imputed data using runMI
outputblimpWLSMV <- cfa.mi(
  model = model,
  data = blimpimp,
  std.lv = TRUE,
  estimator = "WLSMV",
  parameterization = "theta"
)

# view results
summary(outputblimpWLSMV)

# get chi-square statistics
lavTestLRT.mi(outputblimpWLSMV,
              test = "D2",
              pool.robust = TRUE,
              asymptotic = TRUE)

# JOC-FIML ----------------------------------------------------------------
# OpenMx
# FIML model specification in OpenMX
mxOption(NULL, "Default optimizer", "SLSQP")  
# CSOLNP has trouble with thresholds
mxOption(NULL, 'Number of Threads', parallel::detectCores()-4)

# residual Variances
resVars <- mxPath(from = c("CGrip", 
                           "LiftD", "GraspD", "PshPlD"),
                  arrows = 2,
                  free = c(rep(TRUE,1),rep(FALSE,3)),
                  values = c(1))

# means
means <- mxPath(from = "one", 
                to = c("CGrip", 
                       "LiftD", "GraspD", "PshPlD", 
                       "PFUNC"),
                arrows = 1,
                free = c(rep(TRUE,1),rep(FALSE,4)),
                values = c(0))

# latent variances
latVars <- mxPath(from = c("PFUNC"),
                  arrows = 2,
                  connect = "unique.pairs",
                  free = FALSE,
                  values = c(1),
                  labels = c("varPF"))

# free all loadings
facLoads <- mxPath(from = "PFUNC",
                   to = c("CGrip", 
                          "LiftD", "GraspD", "PshPlD"), 
                   arrows = 1,
                   free = c(TRUE))

# set thresholds
thresholds <- mxThreshold(vars = c("LiftD", "GraspD", "PshPlD"), 
                          nThresh = c(2), 
                          free = c(TRUE))

# analyzing the simulated dataset with OpenMx
dataRaw <- mxData(mydata, type = "raw")
OMXmodel <- mxModel("Simulated",
                    type = "RAM",
                    manifestVars = c("CGrip", 
                                     "LiftD", "GraspD", "PshPlD"),
                    latentVars = c("PFUNC"),
                    dataRaw, resVars, latVars, facLoads, means, thresholds)
# fitting hypothesized model
factorFit <- mxRun(OMXmodel)
# fitting saturated and independence model
satFit <- omxSaturatedModel(factorFit, run = TRUE)

# Comparing Parameter Estimates -------------------------------------------

###        get FCS-WLSMV parameter estimates
parEst <- summary(outputMICE, standardized = TRUE)
# only keep the estimates and SE
parEst <- parEst[,grep("lhs|op|rhs|est|se", 
                       names(parEst))]
# create new columns for OpenMx estimates
parEst$est.OMX <- 999
parEst$se.OMX <- 999
# keep only rows regarding factor loadings
parEst.FCS.WLSMV <- subset(parEst,
                           parEst$op == "=~")

# get EMB-WLSMV parameter estimates
parEst <- summary(outputAmelia, standardized = TRUE)
# only keep the estimates and SE
parEst <- parEst[,grep("lhs|op|rhs|est|se", 
                       names(parEst))]
# create new columns for OpenMx estimates
parEst$est.OMX <- 999
parEst$se.OMX <- 999
# keep only rows regarding factor loadings
parEst.EMB.WLSMV <- subset(parEst,
                           parEst$op == "=~")

# get Blimp-WLSMV parameter estimates
parEst <- summary(outputblimpWLSMV, standardized = TRUE)
# only keep the estimates and SE
parEst <- parEst[,grep("lhs|op|rhs|est|se", 
                       names(parEst))]
# create new columns for OpenMx estimates
parEst$est.OMX <- 999
parEst$se.OMX <- 999
# keep only rows regarding factor loadings
parEst.outputblimpWLSMV <- subset(parEst,
                                  parEst$op == "=~")

###        save summary of OpenMx as a separate object
sumOMX <- summary(factorFit)

# combine all factor loadings and se into one dataframe
parEst.all <- cbind(parEst.FCS.WLSMV[,1:5],
                    parEst.EMB.WLSMV[,4:5],
                    parEst.outputblimpWLSMV[,4:5])
# rename column names
names(parEst.all) <- c("lhs",
                       "op",
                       "rhs",
                       "est.FCS.WLSMV",
                       "se.FCS.WLSMV",
                       "est.EMB.WLSMV",
                       "se.EMB.WLSMV",
                       "est.FCSLV.WLSMV",
                       "se.FCSLV.WLSMV")

# add JOC-FIML               
parEst.all$est.JOC.FIML[1:4] <- sumOMX[["parameters"]]$Estimate[1:4]
parEst.all$se.JOC.FIML[1:4] <- sumOMX[["parameters"]]$Std.Error[1:4]

write.csv(data.frame(parEst.all),
          "parEst.all.csv")
