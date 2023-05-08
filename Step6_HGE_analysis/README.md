# Step 6: HGE analysis

## Compare oCGI species patterns in HGEs and non-HGEs using a resampling test
HGE_analysis.sh  

Calls the following:  
&emsp;&emsp;makeHGEtable_quant.py - make table with info on each HGE/non-HGE in each tissue and time point  
&emsp;&emsp;HGEs_resampling.R - perform resampling test on the cluster  
&emsp;&emsp;summarizeResamplingHGEs.py - summarize results of resampling test for further use in R
