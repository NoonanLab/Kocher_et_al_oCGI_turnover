# 10/26/22
# Purpose: make job files to run ChIP-seq pipeline for hs754 data

import sys

inSampleInfo = open(sys.argv[1],'rt')

FilePath_Dict = {}
# FilePath_Dict[set1/set2] = path to Fasta
FilePath_Dict['set1'] = '/ycga-gpfs/sequencers/illumina/sequencerC/runs/220215_A01519_0072_BH7HV3DSX3/Data/Intensities/BaseCalls/Unaligned-150/Project_Ak2267/'
FilePath_Dict['set2'] = '/ycga-gpfs/sequencers/illumina/sequencerD/runs/220328_A00124_0413_AHFMJFDSX3/Data/Intensities/BaseCalls/Unaligned-150/Project_Ak2267/'

Set_Dict = {'A':'set1', 'B':'set1', 'C':'set1', 'D':'set2', 'E':'set2'}
# Set_Dict[setName] = set1/set2

scratchPath = '/home/ak2267/scratch60/hs754_ChIP/'
projectPath = '/home/ak2267/project/hs754_ChIP/'

outFastqc = open('221026_fastqc_jobFile.txt', 'wt')
outAlign = open('221026_align_jobFile.txt', 'wt')
outSort = open('221026_sort_jobFile.txt', 'wt')
outFilter = open('221026_filter_jobFile.txt', 'wt')
outCallPeaks = open('221026_callPeaks_jobFile.txt', 'wt')
outVisualize = open('221026_visualize_jobFile.txt', 'wt')
outNameSort = open('221026_nameSort_jobFile.txt', 'wt')
outCount = open('221026_count_jobFile.txt', 'wt')

# store info from sampleInfo file to translate YCGA sample names into intuitive names for my use
Sample_Dict = {}
# Sample_Dict[sampleName] = [time point, sampleType, genotype, replicate, filePath]

for line in inSampleInfo:
    splitLine = line.strip().split()
    sampleName = splitLine[0]
    timePoint = splitLine[1]
    sampleType = splitLine[2]
    genotype = splitLine[3]
    replicate = splitLine[4]
    multiplexGroup = sampleName[4:5]
    filePath = FilePath_Dict[Set_Dict[multiplexGroup]]
    Sample_Dict[sampleName] = [timePoint, sampleType, genotype, replicate, multiplexGroup, filePath]
inSampleInfo.close()

# go through Sample_Dict and print to each output file for each sample
# change names to intuitive strings: timePoint_sampleType_genotype_replicate

for sample in Sample_Dict:
    #print(sample)

    timePoint = Sample_Dict[sample][0]
    sampleType = Sample_Dict[sample][1]
    if sampleType == 'Input':
        sampleType = 'input'
    genotype = Sample_Dict[sample][2]
    replicate = Sample_Dict[sample][3]
    multiplexGroup = Sample_Dict[sample][4]
    filePath = Sample_Dict[sample][5]
    
    newNameString = str(multiplexGroup+'_e'+timePoint+'.5_'+sampleType+'_'+genotype+'_'+replicate)

    # make FASTQC job file
    outFastqc.write('cd '+projectPath+'fastqc; module load FastQC/0.11.9-Java-11; fastqc '+filePath+'Sample_'+sample+'_*/'+sample+'*R1_001.fastq.gz '+filePath+'Sample_'+sample+'_*/'+sample+'*R2_001.fastq.gz -o '+projectPath+'fastqc/'+'\n')

    # make ALIGN job file (bowtie2)
    genotype = Sample_Dict[sample][2]
    if genotype == 'WT':
        indexString = '/home/ak2267/genomes/mm39'
    elif genotype == 'HUM':
        indexString = '/home/ak2267/genomes/mm39_humanizedHs754'
    outAlign.write("cd "+scratchPath+"bam ; source /home/ak2267/.bash_profile ; source /home/ak2267/.bashrc ; module load Bowtie2/2.4.2-GCCcore-10.2.0 ; module load SAMtools/1.13-GCCcore-10.2.0 ; bowtie2 --sensitive --no-unal -p 20 -x "+indexString+" -1 "+filePath+"Sample_"+sample+"_*/"+sample+"_*R1_001.fastq.gz -2 "+filePath+"Sample_"+sample+"_*/"+sample+"_*R2_001.fastq.gz "+"2> log/"+newNameString+"_bowtie2.log | samtools view -@7 -b > "+newNameString+"_bowtie2.bam"+'\n')
    # IF RE-RUN, change p to 8
    
    # make SORT job file (samtools)
    outSort.write("cd "+scratchPath+"bam/ ; source /home/ak2267/.bash_profile ; source /home/ak2267/.bashrc ; module load SAMtools/1.13-GCCcore-10.2.0 ; samtools sort -t 2 -o "+newNameString+"_bowtie2_sorted.bam "+scratchPath+"bam/"+newNameString+"_bowtie2.bam"+'\n')
    
    # make FILTER job file - 1st step = remove multimappers and unmapped; 2nd step = remove duplicates
    outFilter.write("cd "+scratchPath+"bam ; source /home/ak2267/.bash_profile ; source /home/ak2267/.bashrc ; module load SAMtools/1.13-GCCcore-10.2.0 ; module load Sambamba/0.7.1 ; sambamba view -h -t 2 -f bam -F \"[XS] == null and not unmapped\" "+newNameString+"_bowtie2_sorted.bam > "+newNameString+"_bowtie2_intermediateFiltered.bam 2> log/"+newNameString+"_intermediateFiltered.log ; sambamba markdup -r -t 2 "+newNameString+"_bowtie2_intermediateFiltered.bam "+newNameString+"_bowtie2_filtered.bam 2> log/"+newNameString+"_filtered.log "+'\n')

    # make CALL PEAKS job file (macs2)
    # H3K27ac/H3K4me3/CTCF are narrow, H3K27me3 is broad (https://www.encodeproject.org/chip-seq/histone/#histone)
    mark = sampleType
    if mark != 'input':
        inputString = multiplexGroup+'_e'+timePoint+'.5_input_'+genotype+'_'+replicate
        
        # narrow for H3K27ac, H3K4me3, CTCF - just remove "--broad", there's no "--narrow" flag
        if mark == 'H3K27ac' or mark == 'H3K4me3' or mark == 'CTCF':
            outCallPeaks.write("cd "+scratchPath+"peaks ; source /home/ak2267/.bash_profile ; source /home/ak2267/.bashrc ; module load MACS2/2.2.7.1-foss-2020b-Python-3.8.6 ; module load SAMtools/1.13-GCCcore-10.2.0 ; macs2 callpeak -t "+scratchPath+"bam/"+newNameString+"_bowtie2_filtered.bam -c "+scratchPath+"bam/"+inputString+"_bowtie2_filtered.bam -n "+newNameString+" -g mm "+'\n')
        # narrow for H3K27me3 - include "--broad"
        if mark == 'H3K27me3':
            outCallPeaks.write("cd "+scratchPath+"peaks ; source /home/ak2267/.bash_profile ; source /home/ak2267/.bashrc ; module load MACS2/2.2.7.1-foss-2020b-Python-3.8.6 ; module load SAMtools/1.13-GCCcore-10.2.0 ; macs2 callpeak -t "+scratchPath+"bam/"+newNameString+"_bowtie2_filtered.bam -c "+scratchPath+"bam/"+inputString+"_bowtie2_filtered.bam -n "+newNameString+" -g mm --broad"+'\n')

    # make VISUALIZE job file
    # for mark files: index bam, make bigWig for visualizing, and convert peak bed file to bigBed
    if mark != 'input':
        peakType = {'H3K4me3':'narrow', 'H3K27ac':'narrow', 'CTCF':'narrow', 'H3K27me3':'broad'}
        outVisualize.write("cd "+scratchPath+"bam ; source /home/ak2267/.bash_profile ; source /home/ak2267/.bashrc ; module load deepTools/3.5.1-foss-2020b ; samtools index "+scratchPath+"bam/"+newNameString+"_bowtie2_filtered.bam ; cd "+scratchPath+"visualize ; bamCoverage -b "+scratchPath+"bam/"+newNameString+"_bowtie2_filtered.bam -o "+scratchPath+"visualize/"+newNameString+".bw --effectiveGenomeSize 2309746861 --normalizeUsing CPM --binSize 10 --extendReads 300 --centerReads -p 2 2> log/"+newNameString+"_bamCoverage.log ; cut -f 1,2,3,4 "+scratchPath+"peaks/"+newNameString+"_peaks."+peakType[mark]+"Peak > temp_"+newNameString+".bed ; bedToBigBed temp_"+newNameString+".bed /home/ak2267/genomes/chrom.sizes/mm39.chrom.sizes "+newNameString+"_"+peakType[mark]+"Peak.bb"+'\n')
        
    # for input files: index bam, just make bigWig for visualizing, as there is no bed file
    elif mark == 'input':
        outVisualize.write('cd '+scratchPath+"bam ; source /home/ak2267/.bash_profile ; source /home/ak2267/.bashrc ; module load deepTools/3.5.1-foss-2020b ; samtools index "+scratchPath+"bam/"+newNameString+"_bowtie2_filtered.bam ; cd "+scratchPath+"visualize ; bamCoverage -b "+scratchPath+"bam/"+newNameString+"_bowtie2_filtered.bam -o "+scratchPath+"visualize/"+newNameString+".bw --effectiveGenomeSize 2309746861 --normalizeUsing CPM --binSize 10 --extendReads 300 --centerReads -p 2 2> log/"+newNameString+"_bamCoverage.log"+'\n')
        
    # add re-sort step to sort the filtered bam files by name for HTSeq - position causes memory fails
    outNameSort.write('cd '+scratchPath+'bam/ ; source /home/ak2267/.bash_profile ; source /home/ak2267/.bashrc ; module load SAMtools/1.13-GCCcore-10.2.0 ; samtools sort -n -o '+newNameString+'_bowtie2_filtered_nameSorted.bam '+scratchPath+'bam/'+newNameString+'_bowtie2_filtered.bam'+'\n')

    # count reads in peaks with HTSeq - requires first making consensus peaks with liftover scheme described in hs754_ChIP-seq.sh
    # actually want to run this on this inputs too, given potential copy number issues on chr19
    # but run for each mark since each mark has different peaks called
    if mark != 'input':
        inputString = multiplexGroup+'_e'+timePoint+'.5_input_'+genotype+'_'+replicate
        fileNameInput = multiplexGroup+'_e'+timePoint+'.5_input'+mark+'_'+genotype+'_'+replicate
        if genotype == 'WT':
            # run for actual peaks
            outCount.write('cd '+projectPath+'counts ; source /home/ak2267/.bash_profile ; source /home/ak2267/.bashrc ; module load HTSeq ; htseq-count -r name -f bam -s no '+scratchPath+'bam/'+newNameString+'_bowtie2_filtered_nameSorted.bam '+projectPath+'mergedInMm39/'+multiplexGroup+'_e'+timePoint+'.5_'+mark+'_mergedPeaks.gtf > '+newNameString+'.quant'+'\n')
            # run for input files
            outCount.write('cd '+projectPath+'counts ; source /home/ak2267/.bash_profile ; source /home/ak2267/.bashrc ; module load HTSeq ; htseq-count -r name -f bam -s no '+scratchPath+'bam/'+inputString+'_bowtie2_filtered_nameSorted.bam '+projectPath+'mergedInMm39/'+multiplexGroup+'_e'+timePoint+'.5_'+mark+'_mergedPeaks.gtf > '+fileNameInput+'.quant'+'\n')
        if genotype == 'HUM':
            # run for actual peaks
            outCount.write('cd '+projectPath+'counts ; source /home/ak2267/.bash_profile ; source /home/ak2267/.bashrc ; module load HTSeq ; htseq-count -r name -f bam -s no '+scratchPath+'bam/'+newNameString+'_bowtie2_filtered_nameSorted.bam '+projectPath+'mergedInHUMmm39/'+multiplexGroup+'_e'+timePoint+'.5_'+mark+'_mergedPeaks_inHUMmm39.gtf > '+newNameString+'.quant'+'\n')
            # run for input files
            outCount.write('cd '+projectPath+'counts ; source /home/ak2267/.bash_profile ; source /home/ak2267/.bashrc ; module load HTSeq ; htseq-count -r name -f bam -s no '+scratchPath+'bam/'+inputString+'_bowtie2_filtered_nameSorted.bam '+projectPath+'mergedInHUMmm39/'+multiplexGroup+'_e'+timePoint+'.5_'+mark+'_mergedPeaks_inHUMmm39.gtf > '+fileNameInput+'.quant'+'\n')
            
outFastqc.close()
outAlign.close()
outSort.close()
outFilter.close()
outCallPeaks.close()
outVisualize.close()
outNameSort.close()
outCount.close()

# MACS2 parameters:
# -t = treatment file
# -c = control file (i.e. input)
# -n = stem of name for output
# --narrow = narrow peaks, use for H3K27ac and H3K4me3 - based on ENCODE standards
# --broad = broad peaks, use for H3K4me1
# -g = genome size, hs indicates to use 2.7e9 - for rest of species use <> to calculate
# --bdg = make bedgraph files as output
