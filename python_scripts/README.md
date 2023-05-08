# This directory contains python scripts that are reused throughout the entire pipeline.

## To restrict bed files to sites that map uniquely between two species:
restrictToSitesThatMapBack.py  
restrictToLO.py

## To sort the names in column 4 of a bed file, given that format involves mix of a_#### and b_####
sortRegionNames.py 

## To make a GTF file from a BED file, taking the "exon" name from column 4 of the bed file - often is a bed file of ChIP peaks, not exons
makeGTF_nameFromCol4.py

## To output average RPM and RPKM across replicates from the output of featureCounts
countRPM_general.py - for all peaks in a file (for use in singleSpecies pipelines)  
countRPM.py - restricts output to reconciledPeaks with names containing a_#### or b_#### (for use in speciesPairs pipelines)

## To output average bigwig signal across replicates from the output of bigWigAverageOverbed. "RPM" and "RPKM" are actually total and average bigWig signal, not literal RPM and RPKM
countBigWigSignal_general.py - for all peaks in a file (for use in singleSpecies pipelines)  
countBigWigSignal.py - restricts output to reconciledPeaks with names containing a_#### or b_#### (for use in speciesPairs pipelines)
