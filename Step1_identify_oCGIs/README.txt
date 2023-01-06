# This directory contains scripts for calling CGIs in each genome
# And for generating feature annotations in each genome, which will be used to filter CGIs to oCGIs
# Genome versions compatible with data from Roller & Stamper et al 2021: rheMac10, calJac4, mm39, rn7, susScr11, canFam6, felCat9, equCab3
# Genome versions comparible with data from Reilly, Yin et al 2015 and Cotney, Leng et al 2013: hg19, rheMac2, mm9

Call CGIs in each genome:
220717_UCSC_AL_CpGislands.sh

Generate bed files for feature annotations in each species (based on RefSeq, FANTOM TSSs, and ENCODE blacklist regions):
220724_annotateFeatures.sh
^Calls the script python_scripts/makeCategoryBedFiles_plusIntrons_chrSize.py - to generate bed files with RefSeq annotated features

220725_rmsk_blacklist.sh
