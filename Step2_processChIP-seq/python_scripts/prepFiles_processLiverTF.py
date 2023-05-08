# 7/20/22
# Purpose: make job files to download & process Liver TF data from Ballester 2014 eLife
# Info from this file https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-1509/E-MTAB-1509.sdrf.txt

# cd /home/ak2267/project/LiverTF/batch
# python prepFiles/prepFiles_processLiverTF.py E-MTAB-1509.sdrf.txt

import sys

inFileSummary = open(sys.argv[1],'rt')

Sample_Dict = {}
# Sample_Dict[species][mark][individual][a/b/c for split files] = fastq

speciesList = ['Mus musculus','Rattus norvegicus','Canis familiaris','Macaca mulatta']
speciesListCommon = ['mouse','rat','dog','macaque']
Species_Dict = {'Mus musculus':'mouse','Rattus norvegicus':'rat','Canis familiaris':'dog','Macaca mulatta':'macaque','Homo sapiens':'human'}

# Build dictionaries with all info

for line in inFileSummary:
    splitLine = line.strip().split('\t')
    if splitLine[0] != 'Source Name':

        #for i in range(0,len(splitLine)):
        #    print(str(i)+'\t'+str(splitLine[i]))
        
        species = str(splitLine[3]).strip()
        speciesCommon = Species_Dict[species]
        tissue = str(splitLine[6])
        individual = str(splitLine[5])

        if str(splitLine[29])[0] == 'i':
            mark = 'input'
        else:
            mark = str(splitLine[29])

        # Sample_Dict[species][mark][individual][a/b/c for split files] = fastq

        if speciesCommon in speciesListCommon:
            fastq = str(splitLine[27])

            if species not in Sample_Dict:
                Sample_Dict[species] = {}
            if mark not in Sample_Dict[species]:
                Sample_Dict[species][mark] = {}
            if individual not in Sample_Dict[species][mark]:
                Sample_Dict[species][mark][individual] = {}
            for letter in ['a','b','c','d','e']:
                if letter not in Sample_Dict[species][mark][individual]:
                    Sample_Dict[species][mark][individual][letter] = fastq
                    break

#for species in Sample_Dict:
 #   for mark in Sample_Dict[species]:
  #      for individual in Sample_Dict[species][mark]:
   #         for letter in Sample_Dict[species][mark][individual]:
                #print(species+'\t'+mark+'\t'+individual+'\t'+letter)

Species_Translator = {'Canis familiaris':'canFam6','Macaca mulatta':'rheMac10','Mus musculus':'mm39','Rattus norvegicus':'rn7'}
Effective_GenomeSize = {'canFam6':2312743346,'rheMac10':2936892725,'mm39':2654621783,'rn7':2626580772}

outDownload = open('220720_downloadFastq_jobFile.txt','wt')
outFastqc = open('220720_fastqc_jobFile.txt','wt')
outAlign = open('220720_align_jobFile.txt','wt')
outSort = open('220720_sort_jobFile.txt','wt')

outCombine = open('220720_combine_jobFile.txt','wt')
outFilter = open('220720_filter_jobFile.txt','wt')
outRename = open('220721_rename_jobFile.txt','wt')
outCallPeaks = open('220721_callPeaks_jobFile.txt','wt')
outVisualize = open('220721_visualize_jobFile.txt','wt')

# Dictionary to rename individuals as replicates (for ease of reading file names)
# Mouse is the tricky species where input replicates don't exist for all individuals
# Split these samples such that there is a rep 1 and rep 2 for each mark based on the individuals, and then
# use mmuBL60ON562 as input 1 and mmuOON489 as input 2
Replicate_Dict = {'cfa3':'1','cfa4':'2','mmlBlues':'1','mmlBob':'2','rno5':'1','rno8':'2','mmu0h490+491':'1','mmu12':'2','mmuBL60ON562':'1','mmuOON404':'1','mmuOON405':'2','mmuOON489':'2'}
# Replicate_Dict[individual] = replicate

for species in Sample_Dict:
    for mark in Sample_Dict[species]:
        for individual in Sample_Dict[species][mark]:
            fileStem = str(Species_Translator[species])+'_'+str(mark)+'_'+str(individual)

            for splitFile in Sample_Dict[species][mark][individual]:
                fileString = fileStem +'_'+ str(splitFile)
                #print(fileString)

                # make DOWNLOAD file for use with deadSimpleQueue (dSQ)
                outDownload.write("cd /gpfs/gibbs/pi/noonan/ak2267/LiverTF/fastq ; source /home/ak2267/.bash_profile ; source /home/ak2267/.bashrc ; wget "+str(Sample_Dict[species][mark][individual][splitFile])+" ; mv "+str(Sample_Dict[species][mark][individual][splitFile].split('/')[-1])+" "+fileString+".fastq.gz"+'\n')

                # make FASTQC batch files
                outFastqc.write("cd /gpfs/gibbs/pi/noonan/ak2267/LiverTF/fastq ; source /home/ak2267/.bash_profile ; source /home/ak2267/.bashrc ; module load FastQC/0.11.9-Java-11 ; fastqc "+fileString+".fastq.gz -o fastqc"+'\n')

                # make ALIGN batch files - remove '--no-unal' and change to '--very-sensitive'
                outAlign.write("cd /gpfs/gibbs/pi/noonan/ak2267/LiverTF/bam ; source /home/ak2267/.bash_profile ; source /home/ak2267/.bashrc ; module load Bowtie2/2.4.2-GCCcore-10.2.0 ; module load SAMtools/1.13-GCCcore-10.2.0 ; bowtie2 --very-sensitive -p 8 -x /gpfs/gibbs/pi/noonan/ak2267/genomes/"+str(Species_Translator[species])+" -U /gpfs/gibbs/pi/noonan/ak2267/LiverTF/fastq/"+fileString+".fastq.gz 2> log/"+fileString+"_bowtie2.log | samtools view -@7 -b > "+fileString+"_bowtie2.bam"+'\n')

                # make SORT batch files to sort bams (then remove bam files to free space)
                outSort.write("cd /gpfs/gibbs/pi/noonan/ak2267/LiverTF/bam ; source /home/ak2267/.bash_profile ; source /home/ak2267/.bashrc ; module load SAMtools/1.13-GCCcore-10.2.0 ; module load Sambamba/0.7.1 ; sambamba sort -t 2 -o "+fileString+"_sorted.bam "+fileString+"_bowtie2.bam"+'\n')

            ## following commands happen on combined files for each replicate
            # combine split files
            splitFileArray = []
            for splitFile in Sample_Dict[species][mark][individual]:
                splitFileString = fileStem+'_'+splitFile+'_sorted.bam'
                splitFileArray.append(splitFileString)
            if len(splitFileArray) == 1:
                outCombine.write("cd /gpfs/gibbs/pi/noonan/ak2267/LiverTF/bam ; source /home/ak2267/.bash_profile ; source /home/ak2267/.bashrc ; mv "+splitFileArray[0]+' '+fileStem+"_sorted.bam "+'\n')
            else:
                outCombine.write("cd /gpfs/gibbs/pi/noonan/ak2267/LiverTF/bam ; source /home/ak2267/.bash_profile ; source /home/ak2267/.bashrc ; module load SAMtools/1.13-GCCcore-10.2.0 ; samtools merge "+fileStem+"_sorted.bam "+' '.join(splitFileArray)+'\n')

            # make FILTER batch files - 1st step = remove multimappers and unmapped; 2nd step = remove duplicates
            outFilter.write("cd /gpfs/gibbs/pi/noonan/ak2267/LiverTF/bam ; source /home/ak2267/.bash_profile ; source /home/ak2267/.bashrc ; module load SAMtools/1.13-GCCcore-10.2.0 ; module load Sambamba/0.7.1 ; sambamba view -h -t 2 -f bam -F \"[XS] == null and not unmapped\" "+fileStem+"_sorted.bam > "+fileStem+"_intermediateFiltered.bam 2> log/"+fileStem+"_intermediateFiltered.log ; sambamba markdup -r -t 2 "+fileStem+"_intermediateFiltered.bam "+fileStem+"_filtered.bam 2> log/"+fileStem+"_filtered.log "+'\n')

            # output lines of code to rename files based on replicates
            fileStemRenamed = Species_Translator[species]+'_'+mark+'_'+Replicate_Dict[individual]
            outRename.write('mv '+fileStem+'_filtered.bam '+fileStemRenamed+'_filtered.bam'+'\n')

            # make CALL PEAKS batch files - includes renaming individuals as replicates
            if mark != 'input':
                # narrow for all TFs - just remove "--broad", there's no "--narrow" flag
                outCallPeaks.write("cd /gpfs/gibbs/pi/noonan/ak2267/LiverTF/peaks ; source /home/ak2267/.bash_profile ; source /home/ak2267/.bashrc ; module load MACS2/2.2.7.1-foss-2020b-Python-3.8.6 ; module load SAMtools/1.13-GCCcore-10.2.0 ; macs2 callpeak -t /gpfs/gibbs/pi/noonan/ak2267/LiverTF/bam/"+fileStemRenamed+"_filtered.bam -c /gpfs/gibbs/pi/noonan/ak2267/LiverTF/bam/"+str(Species_Translator[species])+'_input_'+str(Replicate_Dict[individual])+"_filtered.bam -n "+fileStemRenamed+" -g "+str(Effective_GenomeSize[Species_Translator[species]])+' 2> log/'+fileStemRenamed+'.log'+'\n')

            # make VISUALIZE job file - makes bigWigs and bigBeds for use on the genome browser
            # added --extendReads and --centerReads, removed --smoothLength 60
            if mark != 'input':
                peakType = {'CEBPA':'narrow','FOXA1':'narrow','HNF4A':'narrow','HNF6':'narrow'}
                outVisualize.write("cd /gpfs/gibbs/pi/noonan/ak2267/LiverTF/bam ; source /home/ak2267/.bash_profile ; source /home/ak2267/.bashrc ; module load SAMtools/1.13-GCCcore-10.2.0 ; module load deepTools/3.5.1-foss-2020b ; samtools index "+fileStemRenamed+"_filtered.bam ; cd /gpfs/gibbs/pi/noonan/ak2267/LiverTF/visualize ; bamCoverage -b /gpfs/gibbs/pi/noonan/ak2267/LiverTF/bam/"+fileStemRenamed+"_filtered.bam -o "+"/gpfs/gibbs/pi/noonan/ak2267/LiverTF/visualize/"+fileStemRenamed+".bw --effectiveGenomeSize "+str(Effective_GenomeSize[Species_Translator[species]])+" --normalizeUsing CPM --binSize 10 --extendReads 300 --centerReads -p 2 2> log/"+fileStemRenamed+"_bamCoverage.log ; cut -f 1,2,3,4 /gpfs/gibbs/pi/noonan/ak2267/LiverTF/peaks/"+fileStemRenamed+"_peaks."+peakType[mark]+"Peak > temp_"+fileStemRenamed+".bed ; bedToBigBed temp_"+fileStemRenamed+".bed /home/ak2267/genomes/chrom.sizes/"+str(Species_Translator[species])+".chrom.sizes "+fileStemRenamed+"_"+peakType[mark]+"Peak.bb"+'\n')
            elif mark == 'input':
                outVisualize.write("cd /gpfs/gibbs/pi/noonan/ak2267/LiverTF/bam ; source /home/ak2267/.bash_profile ; source /home/ak2267/.bashrc ; module load SAMtools/1.13-GCCcore-10.2.0 ; module load deepTools/3.5.1-foss-2020b ; samtools index "+fileStemRenamed+"_filtered.bam ; cd /gpfs/gibbs/pi/noonan/ak2267/LiverTF/visualize ; bamCoverage -b /gpfs/gibbs/pi/noonan/ak2267/LiverTF/bam/"+str(Species_Translator[species])+'_input_'+str(Replicate_Dict[individual])+"_filtered.bam -o "+"/gpfs/gibbs/pi/noonan/ak2267/LiverTF/visualize/"+str(Species_Translator[species])+'_inputTF_'+str(Replicate_Dict[individual])+".bw --effectiveGenomeSize "+str(Effective_GenomeSize[Species_Translator[species]])+" --normalizeUsing CPM --binSize 10 --extendReads 300 --centerReads -p 2 2> log/"+str(Species_Translator[species])+'_input_'+str(Replicate_Dict[individual])+"_bamCoverage.log"+'\n')
            
            
outDownload.close()
outFastqc.close()
outAlign.close()
outSort.close()
outCombine.close()
outFilter.close()
outRename.close()
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


# MAKE trackDb_${species}_LiverTF.txt FILES TO HOST VISUALIZATIONS ON LAB SERVER

## speciesList = ['rheMac10','mm39','rn7','canFam6']
## markList = ['CEBPA','FOXA1','HNF4A','HNF6','input']

## FilesExist_Dict = {}
## for species in speciesList:
##     FilesExist_Dict[species] = {}
## FilesExist_Dict['rheMac10'] = {'CEBPA':1,'FOXA1':1,'HNF4A':1,'HNF6':2,'input':2}
## FilesExist_Dict['mm39'] = {'CEBPA':2,'FOXA1':2,'HNF4A':2,'HNF6':2,'input':2}
## FilesExist_Dict['rn7'] = {'CEBPA':2,'FOXA1':2,'HNF4A':2,'HNF6':2,'input':2}
## FilesExist_Dict['canFam6'] = {'CEBPA':2,'FOXA1':2,'HNF4A':2,'HNF6':1,'input':2}

## for species in speciesList:
##     if species == 'mm39':
##         count = 113
##     else:
##         count = 85
        
##     outTrackDb = open('trackDb_'+species+'_LiverTFs.txt','wt')

##     #outGenomes.write('genome '+species+'\n')
##     #outGenomes.write('trackDb ./'+species+'/trackDb_'+species+'.txt'+'\n'+'\n')

##     # write to trackDb_{species}.txt files

##     # https://www.color-hex.com/color-palette/1294
##     colorDict = {'CEBPA':'1,31,75','FOXA1':'0,91,150','HNF4A':'100,151,177','HNF6':'111,168,220','inputTF':'209,209,209'}
##     peakType = {'CEBPA':'narrow','FOXA1':'narrow','HNF4A':'narrow','HNF6':'narrow'}
##     trackHeight = {'CEBPA':'5','FOXA1':'5','HNF4A':'5','HNF6':'5','inputTF':'5'}
    
##     for mark in markList:
        
##         if FilesExist_Dict[species][mark] == 2:
##             replicateList = ['1','2']
##         if FilesExist_Dict[species][mark] == 1:
##             replicateList = ['1']

##         if mark == 'input':
##             mark = 'inputTF'
        
##         for replicate in replicateList:

##             # do bigWigs for every mark including input
##             outTrackDb.write('track liver_'+mark+'_'+replicate+'_signal'+'\n')
##             outTrackDb.write('type bigWig 0 '+trackHeight[mark]+'\n')
##             outTrackDb.write('shortLabel liver_'+mark+'_'+replicate+'\n')
##             outTrackDb.write('longLabel liver_'+mark+'_'+replicate+'\n')
##             outTrackDb.write('bigDataUrl '+species+'_'+mark+'_'+replicate+'.bw'+'\n')
##             outTrackDb.write('visibility full'+'\n')
##             outTrackDb.write('color '+colorDict[mark]+'\n')
##             outTrackDb.write('priority '+str(count)+'\n'+'\n')

##             count += 1

##             # do bigBeds only for non-input marks (input has no peak calls and therefore no bigBeds)
##             if mark != 'inputTF':
##                 outTrackDb.write('track liver_'+mark+'_'+replicate+'_peaks'+'\n')
##                 outTrackDb.write('type bigBed'+'\n')
##                 outTrackDb.write('shortLabel liver_'+mark+'_'+replicate+'\n')
##                 outTrackDb.write('longLabel liver_'+mark+'_'+replicate+'\n')
##                 outTrackDb.write('bigDataUrl '+species+'_'+mark+'_'+replicate+'_'+peakType[mark]+'Peak.bb'+'\n')
##                 outTrackDb.write('visibility dense'+'\n')
##                 outTrackDb.write('color '+colorDict[mark]+'\n')
##                 outTrackDb.write('priority '+str(count)+'\n'+'\n')

##                 count += 1

##     outTrackDb.close()

