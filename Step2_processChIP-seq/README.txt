# Step 2: process published ChIP-seq data
Histone modification data from Roller, Stamper et al 2021:  
&emsp;&emsp;map to genomes, remove duplicates, call peaks, identify reproducible peaks, make files for visualization
Histone modiication data from Reilly, Yin et al 2015 and Cotney, Leng et al 2013:  
&emsp;&emsp;download bed files, download bigWig files, identify reproducible peaks

## Process data from Roller, Stamper et al 2021 (adult tissue from brain/liver/muscle/testis, ChIP-seq for H3K4me3, H3K27ac, H3K4me1)
processRollerData_finalRun.sh  

Calls the following python scripts that make slurm job files for batch submission and processing:  
&emsp;&emsp;prepFiles_bowtieIndexes_finalRun.py - downloads genomes and makes bowtie indexes
&emsp;&emsp;prepFiles_processRoller_finalRun.py - downloads Fastq files from ArrayExpress, aligns using bowtie2, sorts & filters multi-mapping and duplicate reads, calls peaks with MACS2, makes bigWigs and bigBeds for visualization
&emsp;&emsp;prepFiles_combineInputs_finalRun.py - combines inputs that come from ArrayExpress as two separate files into single files
&emsp;&emsp;makeIntersectCommands.py # makes bed files with the intersection of all replicate peaks for use downstream

## Process data from Reilly, Yin et al 2015 (developing brain, H3K27ac and H3K4me2) and Cotney, Leng et al 2013 (developing limb, H3K27ac)
processNoonan.sh

## Process TF binding data (CTCF ChIP-seq from Schmidt et al 2012 Cell, focus on rhesus/mouse/rat/dog)
CTCF_dataProcessing.sh  
&emsp;&emsp;Calls prepFiles_processCTCF.py

## Process TF binding data (CEBPA, HNF4a, ONECUT1, and FOXA1 from Ballester et al 2014 eLife, focus on rhesus/mouse/rat/dog)
liverTF_dataProcessing.sh
&emsp;&emsp;Calls prepFiles_processLiverTF.py



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
