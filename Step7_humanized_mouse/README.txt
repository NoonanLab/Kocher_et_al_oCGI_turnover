# Files used in analysis of data generated from humanized mouse line experiments


###### ChIP-seq ######

hs754_ChIP-seq.sh
# Analysis of ChIP-seq from hs754

prepFiles_bowtieIndexesHs754.py
# Generate bowtie indexes for WT mouse and humanized mouse (by swapping human sequence into mm39 fasta file)

prepFiles_hs754_ChIP-seq.py
# Make job files to align, call peaks, and make visualization files for ChIP-seq data

makeTrackHubFiles.py
# Make trackDb.txt file for hosting hs754 data on the browser

# Python scripts for the liftOver pipeline involved in aligning WT and humanized peaks ## NEED TO ORGANIZE THIS
makeHumanizedFasta.py
prepFiles_hs754_ChIP-seq.py
makeTrackHubFiles.py
summarizePeakCategories.py
moveMm39_toHg38.py
adjustHUMpeaksInMm39.py
makeGTF_newPeakNumbers.py
moveHg38_toHUMmm39.py
adjustMergedPeaks_inHUMmm39.py
makeGTF_nameFromCol4.py

###### RNA-seq ######

hs754_RNA-seq.sh
# Analysis of RNA-seq from hs754

prepFiles_hs754_RNA-seq.py
# Make job files to align (and count, included in STAR command) RNA-seq data

