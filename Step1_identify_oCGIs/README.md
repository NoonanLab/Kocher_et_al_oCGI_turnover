# Step 1: identify oCGIs
Each .sh file contains commands used for each processing step  
Each .py file is a python script that performs a step in the pipeline (referenced & described in the .sh files)

## Call CGIs in each genome:
UCSC_AL_CpGislands.sh  
&emsp;&emsp;(Calls the script python_scripts/putCpGsInCol4.py which takes the results of faCount to put the number of CpGs in a CGI in column 4 of a bed file)

Genome versions used for integrating with data from Roller & Stamper et al 2021:  
&emsp;&emsp;rheMac10, calJac4, mm39, rn7, susScr11, canFam6, felCat9, equCab3  
Genome versions used for integrating with data from Reilly, Yin et al 2015 and Cotney, Leng et al 2013:  
&emsp;&emsp;hg19, rheMac2, mm9

## Generate bed files with feature annotations (RefSeq for all genomes, FANTOM & ENCODE blacklist for human and mouse) for filtlering CGIs to oCGIs  
RefSeq: annotateRefSeqFeatures.sh  
&emsp;&emsp;(Calls the script python_scripts/makeCategoryBedFiles_plusIntrons_chrSize.py to generate bed files with RefSeq annotated features)
FANTOM: annotateFANTOMfeatures.sh  
ENCODE blacklist: annotateENCODEblacklist.sh

## Download phastCons, age segmentation (Emera, Yin et al 2016), and RepeatMasker data for use in downstream analysis
PhastCons and age segmentation: phastCons_ageSegmentation.sh
RepeatMasker: rmsk.sh
