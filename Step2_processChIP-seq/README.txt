# Purpose of the scripts in this directory:
# Process ChIP-seq data from Roller, Stamper et al 2021 - map to genomes, remove duplicates, call peaks, identify reproducible peaks, make files for visualization
# Process ChIP-seq data from Reilly, Yin et al 2015 and Cotney, Leng et al 2013 - download bed files, identify reproducible peaks

###############
Process data from Roller, Stamper et al 2021 (adult tissue from brain/liver/muscle/testis, ChIP-seq for H3K4me3, H3K27ac, H3K4me1)
220711_processRollerData_finalRun.sh

Calls the following python scripts that make slurm job files for batch submission and processing:
prepFiles_bowtieIndexes_finalRun.py # download genomes and makes bowtie indexes
prepFiles_processRoller_finalRun.py # download Fastq files from ArrayExpress, aligns using bowtie2, sorts & filters multi-mapping and duplicate reads, calls peaks with MACS2, makes bigWigs and bigBeds for visualization
prepFiles_combineInputs_finalRun.py # combine inputs that come from ArrayExpress as two separate files
makeIntersectCommands.py # make bed files with the intersection of all replicate peaks for use downstream

###############
Process data from Reilly, Yin et al 2015 (developing brain, H3K27ac and H3K4me2) and Cotney, Leng et al 2013 (developing limb, H3K27ac)
220628_processNoonan.sh

###############
Process TF binding data
CTCF ChIP-seq from Schmidt et al 2012 Cell
220615_CTCF_dataProcessing.sh
Calls:
prepFiles_processCTCF.py

Ballester et al 2014 
220720_liverTF_dataProcessing.sh
Calls:
prepFiles_processLiverTF.py
