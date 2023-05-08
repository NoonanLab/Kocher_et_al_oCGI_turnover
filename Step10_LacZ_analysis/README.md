# Step 10: LacZ analysis
Determining whether histone modification peaks predict LacZ reporter activity (Fig S2)

## Generate summary tables for use in R:
LacZ_activity.sh  

Calls the following scripts:  
makeManyBedsFromVISTA.py  
&emsp;&emsp;Takes as input the file "2022-07-15_VISTA.tsv"  
&emsp;&emsp;Outputs a bed file with tested LacZ fragments for each combination of time point and vector backbone  
&emsp;&emsp;For my analysis, wanted to use only fragments tested at E11.5 in the original hsp68 vector  
makeSummaryTableForR_VISTA-centric.py  
&emsp;&emsp;Makes a summary table for downstream R analysis for every combination of  
&emsp;&emsp;histone modification and tissue  

And the following files provided by Dr. Len Pennacchio:  
2022-07-15_VISTA.tsv  
&emsp;&emsp;Summary of the current VISTA database  
tissue.dictionary.csv  
&emsp;&emsp;File for interpreting shorthand codes for tissues where enhancers were active  
