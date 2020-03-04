# Results
This folder contains the dataframes used to generate the diagrams in the Integrated Thesis.
The data are stored in RDS files.

 [1] "mean.o1"     
 [2] "mean.o2"              
 [3] "mean.o3"              
 [4] "mean.c1"              
 [5] "mean.c2"              
 [6] "mean.c3"      
 
 The above columns are the mean factor loadings of o1/o2/o3/c1/c2/c3 across 800 replications.
 
 [7] "mean.o1.error"        
 [8] "mean.o2.error"        
 [9] "mean.o3.error"        
[10] "mean.c1.error"        
[11] "mean.c2.error"        
[12] "mean.c3.error"      

The above columns are the mean error variance of o1/o2/o3/c1/c2/c3 across 800 replications.

[13] "mean.o1.SE"           
[14] "mean.o2.SE"           
[15] "mean.o3.SE"           
[16] "mean.c1.SE"           
[17] "mean.c2.SE"           
[18] "mean.c3.SE"           

The above columns are the mean standard errors of the factor loadings of o1/o2/o3/c1/c2/c3 across 800 replications.

[19] "mean.o1.errorSE"      
[20] "mean.o2.errorSE"      
[21] "mean.o3.errorSE"      
[22] "mean.c1.errorSE"      
[23] "mean.c2.errorSE"      
[24] "mean.c3.errorSE"      

The above columns are the mean standard errors of the error variances of o1/o2/o3/c1/c2/c3 across 800 replications.

[25] "mean.Convergence"     

This column contains the mean number of convergence cases (out of 50) across 800 replications.

[26] "sd.o1"                
[27] "sd.o2"                
[28] "sd.o3"                
[29] "sd.c1"                
[30] "sd.c2"                
[31] "sd.c3"                

The above columns are the standard deviation of the factor loadings of o1/o2/o3/c1/c2/c3 across 800 replications.

[32] "sd.o1.error"          
[33] "sd.o2.error"          
[34] "sd.o3.error"          
[35] "sd.c1.error"          
[36] "sd.c2.error"          
[37] "sd.c3.error"          

The above columns are the standard deviation of the error variance of o1/o2/o3/c1/c2/c3 across 800 replications.

[38] "sd.o1.SE"             
[39] "sd.o2.SE"             
[40] "sd.o3.SE"             
[41] "sd.c1.SE"             
[42] "sd.c2.SE"             
[43] "sd.c3.SE"             

The above columns are the standard deviation of the standard errors of the factor loadings of o1/o2/o3/c1/c2/c3 across 800 replications.

[44] "sd.o1.errorSE"        
[45] "sd.o2.errorSE"        
[46] "sd.o3.errorSE"        
[47] "sd.c1.errorSE"        
[48] "sd.c2.errorSE"        
[49] "sd.c3.errorSE"        

The above columns are the standard errors of the standard errors of the error variances of o1/o2/o3/c1/c2/c3 across 800 replications.

[50] "sd.Convergence"      

This column contains the standard deviation of the number of convergence cases (out of 50) across 800 replications.

[51] "numIter"              

Number of replications ran for each condition, should always be 800.

[52] "method"               

Refers the the "approaches" in the simulation study.

[53] "model"                

Model refers to whether the measurement model was mixed (ordinal + continuous) or purely continuous or ordinal.

[54] "missingProp"     

Factor of simulation study: Missing proportion

[55] "numCat"            

Factor of simulation study: Number of categories

[56] "aSym"      

Factor of simulation study: Degree of asymmetry

[57] "sampleSize"           

Factor of simulation study: Sample size

[58] "Nonconvergence.800"   

Number of non-convergence replications out of 800. Includes improper solution as defined in the paper.

[59] "D3chirej"             

Chi-square rejection rates of the *naive* chi-square from the *D3* method.

[60] "noFitD3"              

Number of replication out of 800 that the *lavTestLRT.mi()* was unable to produce model fit statistics.

[61] "D3scaled.chirej"      

Chi-square rejection rates of the *scaled* chi-square from the *D3* method.

[62] "noFitScaleD3"         

Number of replication out of 800 that the *lavTestLRT.mi()* was unable to produce model fit statistics.

[63] "D2chireg"             

Chi-square rejection rates of the *naive* chi-square from the *D2* method.

[64] "noFitD2"              

Number of replication out of 800 that the *lavTestLRT.mi()* was unable to produce model fit statistics.

[65] "D2scaled.chireg"      

Chi-square rejection rates of the *scaled* chi-square from the *D2* method.

[66] "noFitScaleD2"         

Number of replication out of 800 that the *lavTestLRT.mi()* was unable to produce model fit statistics.


[67] "D2scaled.chiregRobust"

Chi-square rejection rates of the *adjusted* chi-square from the *D2* method.

[68] "noFitScaleD2Robust"   

Number of replication out of 800 that the *lavTestLRT.mi()* was unable to produce model fit statistics.
