# 10/26/22
# Purpose: make job files for analyzing hs754 RNA-seq data (e11.5 and e17.5 diencephalon)

import sys

# 10/26/22
# Purpose: make job files to run ChIP-seq pipeline for hs754 data

inSampleInfo = open(sys.argv[1],'rt')

filePath = '/ycga-gpfs/sequencers/illumina/sequencerD/runs/220404_A00124_0419_BHFNHJDSX3/Data/Intensities/BaseCalls/Unaligned-150/Project_Ak2267/'

scratchPath = '/home/ak2267/scratch60/hs754_RNA/'
projectPath = '/home/ak2267/project/hs754_RNA/'

outFastqc = open('221027_fastqc_jobFile.txt', 'wt')
outAlign = open('221027_align_jobFile.txt', 'wt')

# store info from sampleInfo file to translate YCGA sample names into intuitive names for my use
Sample_Dict = {}
# Sample_Dict[sampleName] = [time point, sampleType, genotype, replicate, filePath]

for line in inSampleInfo:
    splitLine = line.strip().split()
    sampleName = splitLine[0]
    timePoint = splitLine[1]
    genotype = splitLine[2]
    replicate = splitLine[3]
    multiplexGroup = sampleName[4:5]
    Sample_Dict[sampleName] = [timePoint, genotype, replicate, multiplexGroup]
inSampleInfo.close()

# go through Sample_Dict and print to each output file for each sample
# change names to intuitive strings: timePoint_sampleType_genotype_replicate

for sample in Sample_Dict:
    #print(sample)

    timePoint = Sample_Dict[sample][0]
    genotype = Sample_Dict[sample][1]
    replicate = Sample_Dict[sample][2]
    multiplexGroup = Sample_Dict[sample][3]
    
    newNameString = str(multiplexGroup+'_e'+timePoint+'.5_'+genotype+'_'+replicate)

    # make FASTQC job file
    outFastqc.write('cd '+projectPath+'fastqc ; module load FastQC/0.11.9-Java-11 ; fastqc '+filePath+'Sample_'+sample+'_*/'+sample+'*R1_001.fastq.gz '+filePath+'Sample_'+sample+'_*/'+sample+'*R2_001.fastq.gz -o '+projectPath+'fastqc/'+'\n')

    # make ALIGN job file (also does COUNT)
    # https://biocorecrg.github.io/RNAseq_course_2019/alnpractical.html
    outAlign.write('cd '+scratchPath+'bam ; source /home/ak2267/.bash_profile ; source /home/ak2267/.bashrc ; module load STAR/2.7.9a-GCCcore-10.2.0 ; STAR --runMode alignReads --runThreadN 8 --readFilesIn '+filePath+'Sample_'+sample+'_*/'+sample+'*R1_001.fastq.gz '+filePath+'Sample_'+sample+'_*/'+sample+'*R2_001.fastq.gz --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate --outFileNamePrefix '+scratchPath+'bam/'+newNameString+' --quantMode GeneCounts --genomeDir '+projectPath+'STARindex/mm39'+'\n')
    

outFastqc.close()
outAlign.close()
