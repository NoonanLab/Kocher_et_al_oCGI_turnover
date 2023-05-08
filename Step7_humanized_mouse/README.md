# Step 7: humanized mouse experiments

## ChIP-seq

### Analyze ChIP-seq data (in directory "ChIP")
hs754_ChIP-seq.sh  

### Calls the following in directory "ChIP/python_scripts":  
&emsp;&emsp;makeHumanizedFasta.py - makes humanized fasta file from mm39  
&emsp;&emsp;prepFiles_bowtieIndexesHs754.py - generates bowtie indexes for WT mouse and humanized mouse  
&emsp;&emsp;prepFiles_hs754_ChIP-seq.py - makes job files to align, call peaks, and make visualization files  
&emsp;&emsp;makeTrackHubFiles.py - makes trackDb.txt file for hosting hs754 data on the browser

### Python scripts used for aligning WT and humanized peaks (in directory "ChIP/python_scripts/alignPeaks"):

summarizePeakCategories.py  
moveMm39_toHg38.py  
adjustHUMpeaksInMm39.py  
makeGTF_newPeakNumbers.py  
moveHg38_toHUMmm39.py  
adjustMergedPeaks_inHUMmm39.py  

### Also uses script in main repository's directory "python_scripts":
makeGTF_nameFromCol4.py

## RNA-seq

### Analyze RNA-seq data (in directory "RNA")
hs754_RNA-seq.sh

### Calls the following in directory "RNA/python_scripts":
&emsp;&emsp;prepFiles_hs754_RNA-seq.py - makes job files to align (and count, included in STAR command)

