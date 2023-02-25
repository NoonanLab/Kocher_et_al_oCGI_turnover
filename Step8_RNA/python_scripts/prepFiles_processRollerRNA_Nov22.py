# 11/21/22
# Purpose: make files for batch submission for downloading, trimming, QC, aligning, and read counting

# cd /home/ak2267/project/Roller/RNA/batch
# python prepFiles_processRollerRNA.py E-MTAB-8122.sdrf.txt

import sys

inFileSummary = open(sys.argv[1],'rt')

Sample_Dict = {}
# Sample_Dict[species][tissue][replicate] = [fastq_1,fastq_2]

for line in inFileSummary:
    if line[0] == 'd':
        splitLine = line.strip().split('\t')
        species = str(splitLine[3])
        tissue = str(splitLine[10].split(" ")[0])
        replicate = str(splitLine[12])
        fastq = str(splitLine[35])

        if species not in Sample_Dict:
            Sample_Dict[species] = {}
        if tissue not in Sample_Dict[species]:
            Sample_Dict[species][tissue] = {}
        if replicate not in Sample_Dict[species][tissue]:
            Sample_Dict[species][tissue][replicate] = [0,0]
        if '_1.fastq.gz' in fastq:
            Sample_Dict[species][tissue][replicate][0] = fastq
        elif '_2.fastq.gz' in fastq:
            Sample_Dict[species][tissue][replicate][1] = fastq

#for i in Sample_Dict:
#    for j in Sample_Dict[i]:
#        for k in Sample_Dict[i][j]:
#            print(str(i)+'\t'+str(j)+'\t'+str(k))        

Species_Translator = {'Callithrix jacchus':'calJac4','Canis lupus familiaris':'canFam6','Equus caballus':'equCab3','Felis catus':'felCat9','Macaca mulatta':'rheMac10','Monodelphis domestica':'monDom5','Mus musculus':'mm39','Oryctolagus cuniculus':'oryCun2','Rattus norvegicus':'rn7','Sus scrofa':'susScr11'}

outDownload = open('221122_downloadFastq_jobFile.txt','wt')
outFastqc = open('221122_fastqc_jobFile.txt','wt')
outAlign = open('221122_align_jobFile.txt','wt') # also does counting

outFilterDups = open('221205_filterDups_jobFile.txt', 'wt')
outSubsample = open('221205_subsample_jobFile.txt', 'wt')
outFilterDups_and_subsample = open('221205_filterDups_and_subsample_jobFile.txt', 'wt')
outFeatureCounts = open('221205_featureCounts_jobFile.txt', 'wt')

for species in Sample_Dict:
    if species != 'Oryctolagus cuniculus' and species != 'Monodelphis domestica' and species != 'Callithrix jacchus' and species != 'Canis lupus familiaris':
        for tissue in Sample_Dict[species]:
            for replicate in Sample_Dict[species][tissue]:

                fileString = str(Species_Translator[species])+'_'+str(tissue)+'_'+str(replicate)
                #print(fileString)

                # make DOWNLOAD file
                outDownload.write("cd /home/ak2267/scratch60/Roller/fastq_RNA ; source /home/ak2267/.bash_profile ; source /home/ak2267/.bashrc ; wget "+str(Sample_Dict[species][tissue][replicate][0])+" ; mv "+str(Sample_Dict[species][tissue][replicate][0].split('/')[-1])+" "+fileString+"_R1.fastq.gz"+'\n')
                outDownload.write("cd /home/ak2267/scratch60/Roller/fastq_RNA ; source /home/ak2267/.bash_profile ; source /home/ak2267/.bashrc ; wget "+str(Sample_Dict[species][tissue][replicate][1])+" ; mv "+str(Sample_Dict[species][tissue][replicate][1].split('/')[-1])+" "+fileString+"_R2.fastq.gz"+'\n')

                # make FASTQC batch file
                outFastqc.write("cd /home/ak2267/scratch60/Roller/fastq_RNA ; source /home/ak2267/.bash_profile ; source /home/ak2267/.bashrc ; module load FastQC/0.11.9-Java-11 ; fastqc "+fileString+"_R1.fastq.gz -o fastqc"+'\n')
                outFastqc.write("cd /home/ak2267/scratch60/Roller/fastq_RNA ; source /home/ak2267/.bash_profile ; source /home/ak2267/.bashrc ; module load FastQC/0.11.9-Java-11 ; fastqc "+fileString+"_R2.fastq.gz -o fastqc"+'\n')

                # make ALIGN batch files
                # use --quantMode geneCounts to also output counts - see STAR manual for details
                # count against ENSEMBL GTFs
                # https://hbctraining.github.io/Intro-to-rnaseq-hpc-O2/lessons/03_alignment.html
                # STAR --runMode alignReads --runThreadN 8 --readFilesIn '+filePath+'Sample_'+sample+'_*/'+sample+'*R1_001.fastq.gz '+filePath+'Sample_'+sample+'_*/'+sample+'*R2_001.fastq.gz --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate --outFileNamePrefix '+scratchPath+'bam/'+newNameString+' --quantMode GeneCounts --genomeDir '+projectPath+'STARindex/mm39'+'\n')
                outAlign.write("cd /home/ak2267/scratch60/Roller/bam_RNA ; source /home/ak2267/.bash_profile ; source /home/ak2267/.bashrc ; module load STAR/2.7.9a-GCCcore-10.2.0 ; STAR --runMode alignReads --runThreadN 8 --readFilesIn /home/ak2267/scratch60/Roller/fastq_RNA/"+fileString+"_R1.fastq.gz /home/ak2267/scratch60/Roller/fastq_RNA/"+fileString+"_R2.fastq.gz --readFilesCommand zcat --limitGenomeGenerateRAM 40000000000 --limitBAMsortRAM 40000000000 --outSAMtype BAM SortedByCoordinate --outFileNamePrefix /home/ak2267/scratch60/Roller/bam_RNA/"+fileString+" --sjdbGTFfile /home/ak2267/project/Roller/RNA/221121_RNA/GTF/"+str(Species_Translator[species])+".ENSEMBLinUCSC.gtf --quantMode GeneCounts --genomeDir /home/ak2267/scratch60/Roller/genomes/STARindexes/"+str(Species_Translator[species])+'\n')

                # make various filtering / subsampling files
                outFilterDups.write("cd /home/ak2267/scratch60/Roller/bam_RNA ; source /home/ak2267/.bash_profile ; source /home/ak2267/.bashrc ; module load SAMtools/1.13-GCCcore-10.2.0 ; module load Sambamba/0.7.1 ; sambamba view -h -t 2 -f bam -F \"[XS] == null and not unmapped\" "+fileString+"Aligned.sortedByCoord.out.bam > "+fileString+"_intermediateFiltered.bam 2> log/"+fileString+"_intermediateFiltered.log ; sambamba markdup -r -t 2 "+fileString+"_intermediateFiltered.bam "+fileString+"_filtered.bam 2> log/"+fileString+"_filtered.log "+'\n')
                
                # use picard to mark duplicates - sambamba manual says it uses the same tools as picard but with no details
                #outFilterDups.write()
                
                outSubsample.write("cd /home/ak2267/scratch60/Roller/bam_RNA ; source /home/ak2267/.bash_profile ; source /home/ak2267/.bashrc ; module load SAMtools/1.13-GCCcore-10.2.0 ; number=$(samtools view -c "+fileString+"_intermediateFiltered.bam) ; fraction=$(echo \"scale=4;20000000 / $number\" | bc) ; if (( $(echo \"$fraction > 0.99999\" | bc -l) )); then fraction=0.99999; fi ; samtools view -bs ${fraction} "+fileString+"_intermediateFiltered.bam > "+fileString+"_subsampledTo20.bam"+'\n')
                
                outFilterDups_and_subsample.write("cd /home/ak2267/scratch60/Roller/bam_RNA ; source /home/ak2267/.bash_profile ; source /home/ak2267/.bashrc ; module load SAMtools/1.13-GCCcore-10.2.0 ; number=$(samtools view -c "+fileString+"_filtered.bam) ; fraction=$(echo \"scale=4;20000000 / $number\" | bc) ; if (( $(echo \"$fraction > 0.99999\" | bc -l) )); then fraction=0.99999; fi ; samtools view -bs ${fraction} "+fileString+"_filtered.bam > "+fileString+"_filtered_subsampledTo20.bam"+'\n')
                
                # make featureCounts file - count each of the three files above ^^ against ENSEMBL GTF files
                outFeatureCounts.write("cd /home/ak2267/scratch60/Roller/bam_RNA ; source /home/ak2267/.bash_profile ; source /home/ak2267/.bashrc ; module load Subread/2.0.0-GCC-7.3.0-2.30 ; featureCounts -s 2 -p -B -t exon -g gene_id -a /home/ak2267/project/Roller/RNA/221121_RNA/GTF/"+str(Species_Translator[species])+".ENSEMBLinUCSC.gtf -o featureCounts/"+fileString+".filtered.counts "+fileString+"_filtered.bam ; featureCounts -s 2 -p -B -t exon -g gene_id -a /home/ak2267/project/Roller/RNA/221121_RNA/GTF/"+str(Species_Translator[species])+".ENSEMBLinUCSC.gtf -o featureCounts/"+fileString+".subsampledTo20.counts "+fileString+"_subsampledTo20.bam ; featureCounts -s 2 -p -B -t exon -g gene_id -a /home/ak2267/project/Roller/RNA/221121_RNA/GTF/"+str(Species_Translator[species])+".ENSEMBLinUCSC.gtf -o featureCounts/"+fileString+".filtered_subsampledTo20.counts "+fileString+"_filtered_subsampledTo20.bam" +'\n')

outDownload.close()
outFastqc.close()
outAlign.close()

outFilterDups.close()
outSubsample.close()
outFilterDups_and_subsample.close()
outFeatureCounts.close()


# STAR options not used (but could be)
# --outFilterMultimapNmax 1
# --outSAMattrRGline ID:CKp25  # apparently required for featureCounts, but I'm having STAR output counts so this isn't necessary
# --outFilterMismatchNmax 2
