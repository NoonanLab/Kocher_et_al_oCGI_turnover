# 6/15/22
# Purpose: make job files to download & process CTCF data from Schmidt et al 2012 Cell

# cd /home/ak2267/project/Roller/ChIP/batch/CTCF
# python prepFiles/prepFiles_processRoller.py E-MTAB-437.sdrf.txt

import sys

inFileSummary = open(sys.argv[1],'rt')

Sample_Dict = {}
# Sample_Dict[species][tissue][mark][replicate][splitFile] = fastq

speciesList = ['Mus musculus','Rattus norvegicus','Canis familiaris','Macaca mulatta']
speciesListCommon = ['mouse','rat','dog','macaque']

for line in inFileSummary:
    splitLine = line.strip().split('\t')
    if splitLine[0] != 'Source Name':

        #for i in range(0,len(splitLine)):
        #    print(str(i)+'\t'+str(splitLine[i]))
        
        species = str(splitLine[2])
        tissue = str(splitLine[4])

        #print(splitLine[29])
        if str(splitLine[29])[0] == 'i':
            mark = 'input'
            #print(mark)
        else:
            mark = str(splitLine[29])

        # this makes sense for all but mouse
        individual = str(splitLine[22])

        description = splitLine[0].split(' ')
        
        if description[0] in speciesListCommon and len(description) == 3 and species in speciesList:
            fastq = str(splitLine[25])
            replicate = str(splitLine[0].split()[2])

            if replicate == '1' or replicate == '2':

                if species not in Sample_Dict:
                    Sample_Dict[species] = {}
                if tissue not in Sample_Dict[species]:
                    Sample_Dict[species][tissue] = {}
                if mark not in Sample_Dict[species][tissue]:
                    Sample_Dict[species][tissue][mark] = {}
                if replicate not in Sample_Dict[species][tissue][mark]:
                    Sample_Dict[species][tissue][mark][replicate] = {}
                if 'a' in Sample_Dict[species][tissue][mark][replicate]:
                    if 'b' in Sample_Dict[species][tissue][mark][replicate]:
                        Sample_Dict[species][tissue][mark][replicate]['c'] = fastq
                    else:
                        Sample_Dict[species][tissue][mark][replicate]['b'] = fastq
                else:
                    Sample_Dict[species][tissue][mark][replicate]['a'] = fastq

for species in Sample_Dict:
    for tissue in Sample_Dict[species]:
        for mark in Sample_Dict[species][tissue]:
            for replicate in Sample_Dict[species][tissue][mark]:
                for splitFile in Sample_Dict[species][tissue][mark][replicate]:
                    print(species+'\t'+tissue+'\t'+mark+'\t'+replicate+'\t'+splitFile)

Species_Translator = {'Canis familiaris':'canFam6','Macaca mulatta':'rheMac10','Mus musculus':'mm39','Rattus norvegicus':'rn7'}
Effective_GenomeSize = {'canFam6':2312743346,'rheMac10':2936892725,'mm39':2654621783,'rn7':2626580772}

outDownload = open('220620_downloadFastq_jobFile.txt','wt')
outFastqc = open('220620_fastqc_jobFile.txt','wt')
outAlign = open('220620_align_jobFile.txt','wt')
outSort = open('220620_sort_jobFile.txt','wt')

outCombine = open('220620_combine_jobFile.txt','wt')
outFilter = open('220620_filter_jobFile.txt','wt')
outCallPeaks = open('220620_callPeaks_jobFile.txt','wt')
outVisualize = open('220724_visualize_jobFile.txt','wt')

for species in Sample_Dict:
    for tissue in Sample_Dict[species]:
        for mark in Sample_Dict[species][tissue]:
            for replicate in Sample_Dict[species][tissue][mark]:
                fileStem = str(Species_Translator[species])+'_'+str(tissue)+'_'+str(mark)+'_'+str(replicate)
                
                for splitFile in Sample_Dict[species][tissue][mark][replicate]:
                    fileString = fileStem +'_'+ str(splitFile)
                    #print(fileString)

                    # make DOWNLOAD file for use with deadSimpleQueue (dSQ)
                    outDownload.write("cd /home/ak2267/scratch60/CTCF/fastq ; source /home/ak2267/.bash_profile ; source /home/ak2267/.bashrc ; wget "+str(Sample_Dict[species][tissue][mark][replicate][splitFile])+" ; mv "+str(Sample_Dict[species][tissue][mark][replicate][splitFile].split('/')[-1])+" "+fileString+".fastq.gz"+'\n')

                    # make FASTQC batch files
                    outFastqc.write("cd /home/ak2267/scratch60/CTCF/fastq ; source /home/ak2267/.bash_profile ; source /home/ak2267/.bashrc ; module load FastQC/0.11.9-Java-11 ; fastqc "+fileString+".fastq.gz -o fastqc"+'\n')

                    # make ALIGN batch files - remove '--no-unal' and change to '--very-sensitive'
                    outAlign.write("cd /home/ak2267/scratch60/CTCF/bam ; source /home/ak2267/.bash_profile ; source /home/ak2267/.bashrc ; module load Bowtie2/2.4.2-GCCcore-10.2.0 ; module load SAMtools/1.13-GCCcore-10.2.0 ; bowtie2 --very-sensitive -p 8 -x /home/ak2267/scratch60/Roller/genomes/"+str(Species_Translator[species])+" -U /home/ak2267/scratch60/CTCF/fastq/"+fileString+".fastq.gz 2> log/"+fileString+"_bowtie2.log | samtools view -@7 -b > "+fileString+"_bowtie2.bam"+'\n')

                    # make SORT batch files to sort bams (then remove bam files to free space)
                    outSort.write("cd /home/ak2267/scratch60/CTCF/bam ; source /home/ak2267/.bash_profile ; source /home/ak2267/.bashrc ; module load SAMtools/1.13-GCCcore-10.2.0 ; module load Sambamba/0.7.1 ; sambamba sort -t 2 -o "+fileString+"_sorted.bam "+fileString+"_bowtie2.bam"+'\n')

                ## following commands happen on combined files for each replicate
                # combine split files
                splitFileArray = []
                for splitFile in Sample_Dict[species][tissue][mark][replicate]:
                    splitFileString = fileStem+'_'+splitFile+'_sorted.bam'
                    splitFileArray.append(splitFileString)
                if len(splitFileArray) == 1:
                    outCombine.write("cd /home/ak2267/scratch60/CTCF/bam ; source /home/ak2267/.bash_profile ; source /home/ak2267/.bashrc ; mv "+splitFileArray[0]+' '+fileStem+"_sorted.bam "+'\n')
                else:
                    outCombine.write("cd /home/ak2267/scratch60/CTCF/bam ; source /home/ak2267/.bash_profile ; source /home/ak2267/.bashrc ; module load SAMtools/1.13-GCCcore-10.2.0 ; samtools merge "+fileStem+"_sorted.bam "+' '.join(splitFileArray)+'\n')
                
                # make FILTER batch files - 1st step = remove multimappers and unmapped; 2nd step = remove duplicates
                outFilter.write("cd /home/ak2267/scratch60/CTCF/bam ; source /home/ak2267/.bash_profile ; source /home/ak2267/.bashrc ; module load SAMtools/1.13-GCCcore-10.2.0 ; module load Sambamba/0.7.1 ; sambamba view -h -t 2 -f bam -F \"[XS] == null and not unmapped\" "+fileStem+"_sorted.bam > "+fileStem+"_intermediateFiltered.bam 2> log/"+fileStem+"_intermediateFiltered.log ; sambamba markdup -r -t 2 "+fileStem+"_intermediateFiltered.bam "+fileStem+"_filtered.bam 2> log/"+fileStem+"_filtered.log "+'\n')

                # make CALL PEAKS batch files - remove '--bdg' for now to save space, modify to use filtered bam files
                if mark != 'input':
                    # narrow for CTCF - just remove "--broad", there's no "--narrow" flag
                    if mark == 'CTCF':
                        outCallPeaks.write("cd /home/ak2267/scratch60/CTCF/peaks ; source /home/ak2267/.bash_profile ; source /home/ak2267/.bashrc ; module load MACS2/2.2.7.1-foss-2020b-Python-3.8.6 ; module load SAMtools/1.13-GCCcore-10.2.0 ; macs2 callpeak -t /home/ak2267/scratch60/CTCF/bam/"+fileStem+"_filtered.bam -c /home/ak2267/scratch60/CTCF/bam/"+str(Species_Translator[species])+'_'+str(tissue)+'_input_'+str(replicate)+"_filtered.bam -n "+fileStem+" -g "+str(Effective_GenomeSize[Species_Translator[species]])+' 2> log/'+fileStem+'.log'+'\n')

                # make VISUALIZE batch file
                # fileStem = str(Species_Translator[species])+'_'+str(tissue)+'_'+str(mark)+'_'+str(replicate)
                # /gpfs/gibbs/pi/noonan/ak2267/LiverTF/bam/CTCF
                # /gpfs/gibbs/pi/noonan/ak2267/LiverTF/peaks/CTCF
                if mark != 'input':
                    peakType = {'CTCF':'narrow'}
                    outVisualize.write("cd /gpfs/gibbs/pi/noonan/ak2267/LiverTF/visualize ; source /home/ak2267/.bash_profile ; source /home/ak2267/.bashrc ; module load SAMtools/1.13-GCCcore-10.2.0 ; module load deepTools/3.5.1-foss-2020b ; bamCoverage -b /gpfs/gibbs/pi/noonan/ak2267/LiverTF/bam/CTCF/"+fileStem+"_filtered.bam -o "+"/gpfs/gibbs/pi/noonan/ak2267/LiverTF/visualize/"+fileStem+".bw --effectiveGenomeSize "+str(Effective_GenomeSize[Species_Translator[species]])+" --normalizeUsing CPM --binSize 10 --extendReads 300 --centerReads -p 2 2> log/"+fileStem+"_bamCoverage.log ; cut -f 1,2,3,4 /gpfs/gibbs/pi/noonan/ak2267/LiverTF/peaks/CTCF/"+fileStem+"_peaks."+peakType[mark]+"Peak > temp_"+fileStem+".bed ; bedToBigBed temp_"+fileStem+".bed /home/ak2267/genomes/chrom.sizes/"+str(Species_Translator[species])+".chrom.sizes "+fileStem+"_"+peakType[mark]+"Peak.bb"+'\n')
                elif mark == 'input':
                    outVisualize.write("cd /gpfs/gibbs/pi/noonan/ak2267/LiverTF/visualize ; source /home/ak2267/.bash_profile ; source /home/ak2267/.bashrc ; module load SAMtools/1.13-GCCcore-10.2.0 ; module load deepTools/3.5.1-foss-2020b ; bamCoverage -b /gpfs/gibbs/pi/noonan/ak2267/LiverTF/bam/CTCF/"+fileStem+"_filtered.bam -o "+"/gpfs/gibbs/pi/noonan/ak2267/LiverTF/visualize/"+fileStem+".bw --effectiveGenomeSize "+str(Effective_GenomeSize[Species_Translator[species]])+" --normalizeUsing CPM --binSize 10 --extendReads 300 --centerReads -p 2 2> log/"+fileStem+"_bamCoverage.log"+'\n')

outDownload.close()
outFastqc.close()
outAlign.close()
outSort.close()
outCombine.close()
outFilter.close()
outCallPeaks.close()
outVisualize.close()

# MACS2 parameters:
# -t = treatment file
# -c = control file (i.e. input)
# -n = stem of name for output
# --narrow = narrow peaks, use for H3K27ac and H3K4me3 - based on ENCODE standards, this is the default
# --broad = broad peaks, use for H3K4me1
# -g = genome size, hs indicates to use 2.7e9 - for rest of species use <> to calculate
# --bdg = make bedgraph files as output
