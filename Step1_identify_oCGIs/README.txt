# 1/6/22
# This directory contains scripts for calling oCGIs in each genome
# Genome versions compatible with data from Roller & Stamper et al 2021: rheMac10, calJac4, mm39, rn7, susScr11, canFam6, felCat9, equCab3
# Genome versions comparible with data from Reilly, Yin et al 2015 and Cotney, Leng et al 2013: hg19, rheMac2, mm9

Call CGIs in each genome:
220717_UCSC_AL_CpGislands.sh

Generate bed files for feature annotations in each species (based on RefSeq, FANTOM TSSs, and ENCODE blacklist regions):
220724_annotateFeatures.sh
220725_rmsk_blacklist.sh
