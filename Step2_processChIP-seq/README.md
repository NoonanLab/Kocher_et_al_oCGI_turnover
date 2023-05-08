# Step 2: process published ChIP-seq data
Histone modification data from Roller, Stamper et al 2021 (adult brain/liver/muscle/testis, for H3K4me3/H3K27ac/H3K4me1)  
Histone modification data from Reilly, Yin et al 2015 (developing cortex, H3K27ac/H3K4me2) and Cotney, Leng et al 2013 (developing limb, H3K27ac)  
Transcription factor data from Schmidt et al 2012 (CTCF, focus on rhesus/mouse/rat/dog)  
Transcription factor data from Ballester et al 2014 (CEBPA, HNF4a, ONECUT1, and FOXA1, focus on rhesus/mouse/rat/dog)

## Process data from Roller, Stamper et al 2021
processRollerData.sh  
&emsp;&emsp;map to genomes, remove duplicates, call peaks, identify reproducible peaks, make files for visualization  

Calls the following python scripts that make slurm job files for batch submission and processing:  
&emsp;&emsp;prepFiles_bowtieIndexes_finalRun.py - downloads genomes and makes bowtie indexes  
&emsp;&emsp;prepFiles_processRoller_finalRun.py - downloads Fastq files from ArrayExpress, aligns using bowtie2, sorts & filters multi-mapping and duplicate reads, calls peaks with MACS2, makes bigWigs and bigBeds for visualization  
&emsp;&emsp;prepFiles_combineInputs_finalRun.py - combines inputs that come from ArrayExpress as two separate files into single files  
&emsp;&emsp;makeIntersectCommands.py - makes bed files with the intersection of all replicate peaks for use downstream  

## Process data from Reilly, Yin et al 2015 and Cotney, Leng et al 2013 (both from Noonan lab)
processNoonanData.sh  
&emsp;&emsp;download bed files, download bigWig files, identify reproducible peaks

## Process TF binding data from Schmidt et al 2012
processCTCFdata.sh  
&emsp;&emsp;Calls prepFiles_processCTCF.py

## Process TF binding data from Ballester et al 2014
processLiverTFdata.sh  
&emsp;&emsp;Calls prepFiles_processLiverTF.py
