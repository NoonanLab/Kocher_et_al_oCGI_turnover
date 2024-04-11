# Step 5: analysis of oCGI turnover and its overlap with histone modification peak turnover across species pairs

## oCGI-CENTRIC ANALYSIS in sub-directory CGI_centric/
speciesPairs_oCGI-centric.sh - generates a table for each species pair x tissue x mark x time point for analysis downstream in R  

Calls several python scripts to make slurm job files:  
&emsp;&emsp;prepFiles_speciesPairs_CGIcentric_Roller.py  
&emsp;&emsp;prepFiles_speciesPairs_CGIcentric_LiverTF.py  
&emsp;&emsp;prepFiles_speciesPairs_CGIcentric_Noonan.py

All of these job files create intermediate files, then call a python script that integrates them into a final table:  
&emsp;&emsp;countFinalOverlaps_speciesPairs_CGIcentric.py

### Run permutation test to assess enrichment in each box in the grid using R on the cluster
runPermutationOnHPC_oCGI-centric.sh  
&emsp;&emsp;Calls speciesPairs_permutation_HPC_oCGI-centric.R

### Analyze permutation results in R - convert p values to q values with correction for multiple testing, and report Fold Differences (FD)
qValues_oCGI-centric.R

### CGI-centric analysis but without step that removes oCGIs overlapping promoter peaks - for use in Fig 2
prepFiles_speciesPairs_CGIcentric_Roller_noPromFilter.py  
prepFiles_speciesPairs_CGIcentric_Noonan_noPromFilter.py

## PEAK-CENTRIC ANALYSIS in sub-directory Peak_centric/
speciesPairs_peak-centric.sh - generates a table for each species pair x tissue x mark x time point for analysis downstream in R

Calls several python scripts to make slurm job files:  
&emsp;&emsp;prepFiles_speciesPairs_peakCentric_Roller.py  
&emsp;&emsp;prepFiles_speciesPairs_peakCentric_Noonan.py

### Run permutation test to assess enrichment in each box in the grid using R on the cluster
runPermutationOnHPC_peak-centric.sh  
&emsp;&emsp;Calls speciesPairs_permutation_HPC_peak-centric.R

### Analyze permutation results in R - convert p values to q values with correction for multiple testing, and report Fold Differences (FD)
qValues_peak-centric.R



## Test minimal difference requirement for CpG Obs/Exp ratio (Figure S30) in sub-directory Minimal_diff_req/
minimalDifferenceRequirement.sh  
&emsp;&emsp;
