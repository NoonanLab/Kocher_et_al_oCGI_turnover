# Step 1: identify oCGIs
Each .sh file contains commands used for each processing step
Each .py file is a python script that performs a step in the pipeline

## Call CGIs in each genome:
220717_UCSC_AL_CpGislands.sh

Genome versions used for integrating with data from Roller & Stamper et al 2021:  
&emsp;&emsp;rheMac10, calJac4, mm39, rn7, susScr11, canFam6, felCat9, equCab3  
Genome versions used for integrating with data from Reilly, Yin et al 2015 and Cotney, Leng et al 2013:  
&emsp;&emsp;hg19, rheMac2, mm9

## Generate bed files with feature annotations (RefSeq for all genomes, FANTOM & ENCODE blacklist for human and mouse) for filtlering CGIs to oCGIs  
RefSeq: 220724_annotateFeatures.sh  
&emsp;&emsp;(Calls the script python_scripts/makeCategoryBedFiles_plusIntrons_chrSize.py to generate bed files with RefSeq annotated features)
FANTOM: 220616_FANTOM_data.sh  
ENCODE blacklist: 220725_blacklist.sh

## Download phastCons, age segmentation (Emera, Yin et al 2016), and RepeatMasker data for use in downstream analysis
PhastCons and age segmentation: 220729_phastCons_ageSegmentation.sh
RepeatMasker: 220725_rmsk.sh
