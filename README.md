# FIML_MI_JOC_MISSINGDATA

1. The Empirical Illustration folder contains the files used for the Empirical Illustration.  
Empirical_Illustration.R is the R script containing the codes to replicate the empirical illustration.  
empiricalblimp.imp is the text file containing the syntax for Blimp.

2. Main_Script.R, genData.R, FCS_WLSMV.R, FCS_MLM.R, EMB_WLSMV.R, EMB_MLM.R, JOC_FIML.R, Extract_WLSMV.R, Extract_MLM.R and Extract_JOC_FIML.R are the R scripts containing the codes used to run the simulation.  
These are not the actual files used to run the simulation as the process was more complicated when the simulation was ran on NUS's High Performance Computing Clusters. However, the commands used are identical to those used in the simulation and provide a reference to the settings and specifications used in the simulation.
 
Main_Script.R is the main body while the other R files contain functions with different specifications depending on the approach, the other R files are loaded according to the approach chosen in Main_Script.R using `source()`.

3. The results folder contains the dataframes used to plot the diagrams in the thesis, they are saved in RDS files.
