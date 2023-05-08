# This directory contains scripts for calling CGIs in each genome
# And for generating feature annotations in each genome, which will be used to filter CGIs to oCGIs (RefSeq, FANTOM, ENCODE blacklist)
# Finally download phastCons and RepeatMasker data for use downstream

# Genome versions used for integrating with data from Roller & Stamper et al 2021: rheMac10, calJac4, mm39, rn7, susScr11, canFam6, felCat9, equCab3
# Genome versions used for integrating with data from Reilly, Yin et al 2015 and Cotney, Leng et al 2013: hg19, rheMac2, mm9

#########
Call CGIs in each genome:
220717_UCSC_AL_CpGislands.sh

#########
Generate bed files for feature annotations in each species (based on RefSeq, FANTOM TSSs, and ENCODE blacklist regions)

RefSeq
220724_annotateFeatures.sh
^Calls the script python_scripts/makeCategoryBedFiles_plusIntrons_chrSize.py - to generate bed files with RefSeq annotated features

FANTOM
220616_FANTOM_data.sh

ENCODE blacklist
220725_blacklist.sh

#########
Download phastCons data and age segmentation data (from Emera, Yin et al 2016) for use downstream
220729_phastCons_ageSegmentation.sh

#########
Download RepeatMasker data for use downstream
220725_rmsk.sh
