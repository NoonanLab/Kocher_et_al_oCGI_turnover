# This directory contains python scripts that are reused throughout the entire pipeline.

# For use in restricting bed files to sites that map uniquely between two species
restrictToSitesThatMapBack.py

# Also for use in restricting bed files to sites that map uniquely between two species
restrictToLO.py

# For sorting the names in column 4 of a bed file, given that format involves mix of a_#### and b_####
sortRegionNames.py

# Making a GTF file from a BED file, taking the "exon" name from column 4 of the bed file - often is a bed file of ChIP peaks, not exons
makeGTF_nameFromCol4.py

# 
countRPM_general.py

# 
countBigWigSignal_general.py
