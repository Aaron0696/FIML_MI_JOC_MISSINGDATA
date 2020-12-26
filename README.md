# FIML_MI_JOC_MISSINGDATA

1. The Empirical Illustration folder contains the files used for the Empirical Illustration.  
Empirical_Illustration.R is the R script containing the codes to replicate the empirical illustration.  
empiricalblimp.imp is the text file containing the syntax for imputation by Blimp.

2. The following scripts are used to replicate the simulation.

* Main.R
* Main_FCSLV.R

The above scripts contain the main body of the simulation, the other .R files below are sourced depending on the approach selected in the main script. There is a separate script for the FCSLV approaches as Blimp was used to impute the data, which require the creation of imputation instructions and exported data.

* genData.R

genData.R contains the function for generating datasets.

* FCS_WLSMV.R
* FCS_FIML.R
* FCS_MLM.R (Not used)
* FCSLV_WLSMV
* FCSLV_FIML
* EMB_WLSMV.R
* EMB_FIML.R
* EMB_MLM.R (Not used)
* JOC_FIML.R

The above scripts contain a function analyze.Data(), which is used to analyze the data generated by genData.R.

* Extract_WLSMV.R
* Extract_MLM.R
* Extract_JOC_FIML.R 

The above scripts contain a function extract(), which is used to transform the outputs created by analze.Data into a dataframe where each row corresponds to one of the 72 conditions in the study.

These are not the actual files used to run the simulation as the process was more complicated when the simulation was ran on NUS's High Performance Computing Clusters. However, the commands used are identical to those used in the simulation and provide a reference to the settings and specifications used in the simulation.
 
3. The results folder contains the dataframes used to plot the diagrams in the thesis, they are saved in RDS files. Each of these RDS contain a R dataframe with 72 rows, with each row corresponding to one of the conditions in the study.
