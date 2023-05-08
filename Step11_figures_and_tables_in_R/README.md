# Step 1: identify oCGIs
Each .sh file contains commands used for each processing step  
Each .py file is a python script that performs a step in the pipeline (referenced & described in the .sh files)

## Call CGIs in each genome:
UCSC_AL_CpGislands.sh  
&emsp;&emsp;Calls putCpGsInCol4.py which takes the results of faCount to put the number of CpGs in a CGI in column 4 of a bed file




Fig1_CGIcentric.R
Fig1_peakCentric.R
## Generate Fig 1, Fig SX-Y

Fig2_orthologous_oCGI.R
## Generate Fig 2B, Fig SX, Table SX

Fig2_speciesComparisons.R
## Generate Fig 2C-E, Fig SX

Fig3-4_grids.R 
## Generate heatmap grids in Fig 3 & 4

Fig3_rhesusMouse_example.R
## Generate numbers and histogram for example in Fig S3.1

Fig3_peakStatComparisons.R
## Comparisons of peak strength, phastCons, age, etc
