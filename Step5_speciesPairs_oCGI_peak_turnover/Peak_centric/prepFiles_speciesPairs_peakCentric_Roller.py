# 11/30/22
# Purpose: make jobFile for doing analysis to make full spreadsheet for speciesPairs, peakCentric
# Where each line in the spreadsheet is a PEAK that maps between the species
# This file only does Roller histone mod data - see other scripts for Liver TF data and Noonan histone mod data
# Modified 9/15/22 to incorporate phastBias

# 

import sys

speciesArray = ['rheMac10','calJac4','mm39','rn7','susScr11','canFam6','felCat9','equCab3']
Tissues = ['brain','liver','muscle','testis']
Marks = ['H3K4me3','H3K27ac','H3K4me1']

SpeciesCombos = []

inReps = open(sys.argv[1],'rt')
Rep_Dict = {}
# Rep_Dict[species_tissue_mark] = [1,2,3] # or subset if a rep is unused
for line in inReps:
    splitLine = line.strip().split()
    Rep_Dict[splitLine[0]+'_'+splitLine[1]+'_'+splitLine[2]] = splitLine[3].split(',')

phastBiasCombos = ['mm39_rn7','canFam6_felCat9']
phastBias_SpeciesTranslator = {'mm39':'mm10', 'rn7':'rn6', 'canFam6':'canFam3', 'felCat9':'felCat8'}

filterFlag = str(sys.argv[2]) # repeatFilter vs noRepeatFilter

# store filePaths and fileEnds
workingDir = '/gpfs/gibbs/pi/noonan/ak2267/speciesPairs/peakCentric'
pathToConsensusCGIs = '/gpfs/gibbs/pi/noonan/ak2267/speciesPairs/consensusCGIs/'
# then: {speciesA}_{speciesB}/{speciesA}_{speciesB}_speciesA_reconciledCGI.bed
#    or {speciesA}_{speciesB}/{speciesA}_{speciesB}_speciesB_reconciledCGI.bed
pathToAnnotations = '/home/ak2267/genomes/RefSeq/featureAnnotations/'
# then: ${species}_allFeatures.bed / _1_proteinCoding_promoters_merged.bed _2_proteinCoding_exons_merged.bed _3_lncRNA_merged.bed _4_ncRNA_merged.bed _5_pseudogene_merged.bed _6_unknown_merged.bed
pathToFANTOM = '/home/ak2267/genomes/FANTOM/'
# then: FANTOM_TSS_${species}.bed
pathToBlacklist = '/home/ak2267/genomes/blacklist/'
# then: blacklist_${species}.bed
pathToPeaks = '/gpfs/gibbs/pi/noonan/ak2267/Roller/peaks/intersection/'
# then: ${species}_${tissue}_${mark}_intersection.bed
pathToPhastBias = '/gpfs/gibbs/pi/noonan/ak2267/speciesPairs/SequenceAnalysis/gBGC_analysis/phastBias_bed/'
# then: ${species}_allChr_phastBias_merged.bed - note this is in hg38 coordinates

FilesToExclude_Dict = {}
FilesToExclude_Dict['filePath'] = {'RefSeq':'/home/ak2267/genomes/RefSeq/featureAnnotations/','FANTOM':'/home/ak2267/genomes/FANTOM/FANTOM_TSS_','blacklist':'/home/ak2267/genomes/blacklist/blacklist_'}
FilesToExclude_Dict['fileEnd'] = {'RefSeq':'_allFeatures.bed','FANTOM':'.bed','blacklist':'.bed'}


for j in range(0,len(speciesArray)):
    for k in range(j+1,len(speciesArray)):
        speciesString = str(speciesArray[j]) +'_'+ str(speciesArray[k])
        SpeciesCombos.append(speciesString)

for speciesCombo in SpeciesCombos:
    for tissue in Tissues:
        for mark in Marks:
            #print(speciesCombo+'\t'+tissue+'\t'+mark)
            speciesA = speciesCombo.split('_')[0]
            speciesB = speciesCombo.split('_')[1]
            reconciledCGI_A = pathToConsensusCGIs+speciesA+'_'+speciesB+'/'+speciesA+'_'+speciesB+'_speciesA_reconciledCGI.bed'
            reconciledCGI_B = pathToConsensusCGIs+speciesA+'_'+speciesB+'/'+speciesA+'_'+speciesB+'_speciesB_reconciledCGI.bed'
            reconciledHuman = pathToConsensusCGIs+speciesA+'_'+speciesB+'/'+speciesA+'_'+speciesB+'_hg38_reconciledCGI.bed'
            phastConsHuman = '/home/ak2267/genomes/phastCons/phastConsElements100way_hg38.bed'
            ageSegmentationHuman = '/home/ak2267/project/EnhancerAge/clade_specific_elements_hg38.bed'

            outputString = []

            #########
            # STEP 0: append info on where to work and source files
            label = speciesCombo +'_'+ tissue +'_'+ mark
            outputString.append('cd '+workingDir+' ; source /home/ak2267/.bashrc ; source /home/ak2267/.bash_profile')
            outputString.append('mkdir Roller_'+label+' ; cd Roller_'+label)

            ######### (analogous to steps 1-3 in consensusCGI files)
            # STEP 1: make files with 4 columns and CGI names (as "peak" names), with 'a' or 'b' preceding in each species
            outputString.append('cat '+pathToPeaks+speciesA+'_'+tissue+'_'+mark+'_intersection.bed | awk \'{ print $1\"\\t\"$2\"\\t\"$3\"\\ta_\"NR }\' > speciesA_peaks_4col.bed')
            outputString.append('cat '+pathToPeaks+speciesB+'_'+tissue+'_'+mark+'_intersection.bed | awk \'{ print $1\"\\t\"$2\"\\t\"$3\"\\tb_\"NR }\' > speciesB_peaks_4col.bed')

            #########
            # STEP 2: get 'orphan' PEAKS in each species that liftover to the other species 
            
            speciesWithFANTOMandBlacklist = ['hg19','mm39','mm9']
            
            # restrict to non-feature (aka orphan) CGIs in both species
            if speciesA in speciesWithFANTOMandBlacklist:
                filesToExcludeA = [FilesToExclude_Dict['filePath']['RefSeq']+speciesA+FilesToExclude_Dict['fileEnd']['RefSeq'], FilesToExclude_Dict['filePath']['FANTOM']+speciesA+FilesToExclude_Dict['fileEnd']['FANTOM'], FilesToExclude_Dict['filePath']['blacklist']+speciesA+FilesToExclude_Dict['fileEnd']['blacklist']]
            else:
                filesToExcludeA = [FilesToExclude_Dict['filePath']['RefSeq']+speciesA+FilesToExclude_Dict['fileEnd']['RefSeq']]

            if speciesB in speciesWithFANTOMandBlacklist:
                filesToExcludeB = [FilesToExclude_Dict['filePath']['RefSeq']+speciesB+FilesToExclude_Dict['fileEnd']['RefSeq'], FilesToExclude_Dict['filePath']['FANTOM']+speciesB+FilesToExclude_Dict['fileEnd']['FANTOM'], FilesToExclude_Dict['filePath']['blacklist']+speciesB+FilesToExclude_Dict['fileEnd']['blacklist']]
            else:
                filesToExcludeB = [FilesToExclude_Dict['filePath']['RefSeq']+speciesB+FilesToExclude_Dict['fileEnd']['RefSeq']]

            hg38FilesToExclude = [FilesToExclude_Dict['filePath']['RefSeq']+'hg38'+FilesToExclude_Dict['fileEnd']['RefSeq'], FilesToExclude_Dict['filePath']['FANTOM']+'hg38'+FilesToExclude_Dict['fileEnd']['FANTOM'], FilesToExclude_Dict['filePath']['blacklist']+'hg38'+FilesToExclude_Dict['fileEnd']['blacklist']]
            hg19FilesToExclude = [FilesToExclude_Dict['filePath']['RefSeq']+'hg19'+FilesToExclude_Dict['fileEnd']['RefSeq'], FilesToExclude_Dict['filePath']['FANTOM']+'hg19'+FilesToExclude_Dict['fileEnd']['FANTOM'], FilesToExclude_Dict['filePath']['blacklist']+'hg19'+FilesToExclude_Dict['fileEnd']['blacklist']]

            # lift both species to human

            # remove features in both species
            outputString.append('bedtools intersect -v -a speciesA_peaks_4col.bed -b '+' -b '.join(filesToExcludeA)+' > speciesA_peaks_4col_noFeature.bed')
            outputString.append('bedtools intersect -v -a speciesB_peaks_4col.bed -b '+' -b '.join(filesToExcludeB)+' > speciesB_peaks_4col_noFeature.bed')

            # store chains
            chainAtoHg38 = '/home/ak2267/genomes/chain/'+speciesA+'ToHg38.over.chain.gz'
            chainBtoHg38 = '/home/ak2267/genomes/chain/'+speciesB+'ToHg38.over.chain.gz'
            chainHg38toA = '/home/ak2267/genomes/chain/hg38To'+speciesA[0].upper()+speciesA[1:]+'.over.chain.gz'
            chainHg38toB = '/home/ak2267/genomes/chain/hg38To'+speciesB[0].upper()+speciesB[1:]+'.over.chain.gz'

            # lift to human
            outputString.append('liftOver -minMatch=0.3 speciesA_peaks_4col_noFeature.bed '+chainAtoHg38+' speciesA_peaks_4col_noFeature_LOhg38.bed unMapped')
            outputString.append('liftOver -minMatch=0.3 speciesB_peaks_4col_noFeature.bed '+chainBtoHg38+' speciesB_peaks_4col_noFeature_LOhg38.bed unMapped')

            # lift back and restrict list in human to only those that lift back to the same CGI
            outputString.append('liftOver -minMatch=0.3 speciesA_peaks_4col_noFeature_LOhg38.bed '+chainHg38toA+' speciesA_peaks_4col_noFeature_liftedBackToA.bed unMapped')
            outputString.append('liftOver -minMatch=0.3 speciesB_peaks_4col_noFeature_LOhg38.bed '+chainHg38toB+' speciesB_peaks_4col_noFeature_liftedBackToB.bed unMapped')
            outputString.append('bedtools intersect -wao -a speciesA_peaks_4col_noFeature_liftedBackToA.bed -b speciesA_peaks_4col_noFeature.bed > speciesA_checkMappingBack.txt')
            outputString.append('bedtools intersect -wao -a speciesB_peaks_4col_noFeature_liftedBackToB.bed -b speciesB_peaks_4col_noFeature.bed > speciesB_checkMappingBack.txt')
            outputString.append('python /gpfs/gibbs/pi/noonan/ak2267/speciesPairs/consensusCGIs/Scripts/restrictToSitesThatMapBack.py speciesA_checkMappingBack.txt > speciesA_peaks_4col_noFeature_mapBack_inAcoord.bed')
            outputString.append('python /gpfs/gibbs/pi/noonan/ak2267/speciesPairs/consensusCGIs/Scripts/restrictToSitesThatMapBack.py speciesB_checkMappingBack.txt > speciesB_peaks_4col_noFeature_mapBack_inBcoord.bed')
            outputString.append('python /gpfs/gibbs/pi/noonan/ak2267/speciesPairs/consensusCGIs/Scripts/restrictToLO.py speciesA_peaks_4col_noFeature_LOhg38.bed speciesA_peaks_4col_noFeature_mapBack_inAcoord.bed > speciesA_peaks_4col_noFeature_LOhg38_mapBack.bed')
            outputString.append('python /gpfs/gibbs/pi/noonan/ak2267/speciesPairs/consensusCGIs/Scripts/restrictToLO.py speciesB_peaks_4col_noFeature_LOhg38.bed speciesB_peaks_4col_noFeature_mapBack_inBcoord.bed > speciesB_peaks_4col_noFeature_LOhg38_mapBack.bed')

            # merge, remove features in human, sort file names to be a_# then b_#
            outputString.append('cat speciesA_peaks_4col_noFeature_LOhg38_mapBack.bed speciesB_peaks_4col_noFeature_LOhg38_mapBack.bed > cat.bed')
            outputString.append('sort -k1,1 -k2,2n cat.bed > sort.bed')
            outputString.append('bedtools merge -c 4 -o collapse -i sort.bed > merge_hg38.bed')
            outputString.append('bedtools intersect -v -a merge_hg38.bed -b '+' -b '.join(hg38FilesToExclude)+' > merge_hg38_noFeature.bed')
            outputString.append('python /gpfs/gibbs/pi/noonan/ak2267/speciesPairs/consensusCGIs/Scripts/sortRegionNames.py merge_hg38_noFeature.bed > merge_hg38_noFeature_sortedNames.bed')

			# lift back out to speciesA and speciesB
            outputString.append('liftOver -minMatch=0.3 merge_hg38_noFeature_sortedNames.bed '+chainHg38toA+' mergedPeaks_liftToA.bed unMapped')
            outputString.append('liftOver -minMatch=0.3 merge_hg38_noFeature_sortedNames.bed '+chainHg38toB+' mergedPeaks_liftToB.bed unMapped')

			# lift back to human and restrict list in speciesA and speciesB to only those that lift back to the same CGI in human (analagous to above)
            outputString.append('liftOver -minMatch=0.3 mergedPeaks_liftToA.bed '+chainAtoHg38+' mergedPeaks_liftToA_backToHg38.bed unMapped')
            outputString.append('liftOver -minMatch=0.3 mergedPeaks_liftToB.bed '+chainBtoHg38+' mergedPeaks_liftToB_backToHg38.bed unMapped')
            outputString.append('bedtools intersect -wao -a mergedPeaks_liftToA_backToHg38.bed -b merge_hg38_noFeature_sortedNames.bed > speciesA_checkMappingBackInHg38.txt')
            outputString.append('bedtools intersect -wao -a mergedPeaks_liftToB_backToHg38.bed -b merge_hg38_noFeature_sortedNames.bed > speciesB_checkMappingBackInHg38.txt')
            outputString.append('python /gpfs/gibbs/pi/noonan/ak2267/speciesPairs/consensusCGIs/Scripts/restrictToSitesThatMapBack.py speciesA_checkMappingBackInHg38.txt > speciesA_merged_thatMapBackToHg38_inHg38.bed')
            outputString.append('python /gpfs/gibbs/pi/noonan/ak2267/speciesPairs/consensusCGIs/Scripts/restrictToSitesThatMapBack.py speciesB_checkMappingBackInHg38.txt > speciesB_merged_thatMapBackToHg38_inHg38.bed')
            outputString.append('python /gpfs/gibbs/pi/noonan/ak2267/speciesPairs/consensusCGIs/Scripts/restrictToLO.py mergedPeaks_liftToA.bed speciesA_merged_thatMapBackToHg38_inHg38.bed > merged_speciesA_mapBack.bed')
            outputString.append('python /gpfs/gibbs/pi/noonan/ak2267/speciesPairs/consensusCGIs/Scripts/restrictToLO.py mergedPeaks_liftToB.bed speciesB_merged_thatMapBackToHg38_inHg38.bed > merged_speciesB_mapBack.bed')

            # remove features
            outputString.append('bedtools intersect -v -a merged_speciesA_mapBack.bed -b '+' -b '.join(filesToExcludeA)+' > mergedPeaks_liftToA_noFeature.bed')
            outputString.append('bedtools intersect -v -a merged_speciesB_mapBack.bed -b '+' -b '.join(filesToExcludeB)+' > mergedPeaks_liftToB_noFeature.bed')

            #########
            # STEP 3: restrict all three lists (in A, in B, and in hg38) to those present in both A and B
            outputString.append('python /gpfs/gibbs/pi/noonan/ak2267/speciesPairs/consensusCGIs/Scripts/restrictToLO.py mergedPeaks_liftToA_noFeature.bed mergedPeaks_liftToB_noFeature.bed > speciesA_reconciledPeaks.bed')
            outputString.append('python /gpfs/gibbs/pi/noonan/ak2267/speciesPairs/consensusCGIs/Scripts/restrictToLO.py mergedPeaks_liftToB_noFeature.bed mergedPeaks_liftToA_noFeature.bed > speciesB_reconciledPeaks.bed')
            outputString.append('python /gpfs/gibbs/pi/noonan/ak2267/speciesPairs/consensusCGIs/Scripts/restrictToLO.py merge_hg38_noFeature_sortedNames.bed mergedPeaks_liftToA_noFeature.bed > intermediateHg38_reconciledFile.bed')
            outputString.append('python /gpfs/gibbs/pi/noonan/ak2267/speciesPairs/consensusCGIs/Scripts/restrictToLO.py intermediateHg38_reconciledFile.bed mergedPeaks_liftToB_noFeature.bed > human_reconciledPeaks.bed')

            #########
            # STEP 4: intersect reconciled PEAKS with several files in each species:
            # peaks to get peak info
            # CGIs - both reconciled and not reconciled
            # rmsk
            # phastCons /home/ak2267/genomes/phastCons
            # age segments

            pythonInputFiles = []

            # intersect each reconciled list with CGI files
            outputString.append('bedtools intersect -wao -a speciesA_reconciledPeaks.bed -b /home/ak2267/genomes/CGI/UCSC_AL/'+speciesA+'_CGIsAL.bed > speciesA_reconciledPeaks_intersectSpeciesA_CGI.txt')
            outputString.append('bedtools intersect -wao -a speciesB_reconciledPeaks.bed -b /home/ak2267/genomes/CGI/UCSC_AL/'+speciesB+'_CGIsAL.bed > speciesB_reconciledPeaks_intersectSpeciesB_CGI.txt')

            pythonInputFiles.append('speciesA_reconciledPeaks_intersectSpeciesA_CGI.txt')
            pythonInputFiles.append('speciesB_reconciledPeaks_intersectSpeciesB_CGI.txt')

            # intersect each reconciled list with RECONCILED CGI files
            outputString.append('bedtools intersect -wao -a speciesA_reconciledPeaks.bed -b '+ pathToConsensusCGIs + speciesCombo +'/'+ speciesCombo + '_speciesA_reconciledCGI.bed > speciesA_reconciledPeaks_intersectSpeciesA_reconciledCGI.txt')
            outputString.append('bedtools intersect -wao -a speciesB_reconciledPeaks.bed -b '+ pathToConsensusCGIs + speciesCombo +'/'+ speciesCombo + '_speciesB_reconciledCGI.bed > speciesB_reconciledPeaks_intersectSpeciesB_reconciledCGI.txt')

            pythonInputFiles.append('speciesA_reconciledPeaks_intersectSpeciesA_reconciledCGI.txt')
            pythonInputFiles.append('speciesB_reconciledPeaks_intersectSpeciesB_reconciledCGI.txt')

            # intersect each reconciled list with RepeatMasker
            outputString.append('bedtools intersect -wao -a speciesA_reconciledPeaks.bed -b /home/ak2267/genomes/rmsk/rmsk_'+speciesA+'_merged.bed > speciesA_reconciledPeaks_intersectSpeciesA_rmsk.txt')
            outputString.append('bedtools intersect -wao -a speciesB_reconciledPeaks.bed -b /home/ak2267/genomes/rmsk/rmsk_'+speciesB+'_merged.bed > speciesB_reconciledPeaks_intersectSpeciesB_rmsk.txt')

            pythonInputFiles.append('speciesA_reconciledPeaks_intersectSpeciesA_rmsk.txt')
            pythonInputFiles.append('speciesB_reconciledPeaks_intersectSpeciesB_rmsk.txt')

            # intersect human reconciled list with phastCons elements
            outputString.append('bedtools intersect -wao -a human_reconciledPeaks.bed -b '+phastConsHuman+' > human_reconciledPeaks_intersectPhastCons.txt')
            pythonInputFiles.append('human_reconciledPeaks_intersectPhastCons.txt')

            # intersect human reconciled list with age segmentation map
            outputString.append('bedtools intersect -wao -a human_reconciledPeaks.bed -b '+ageSegmentationHuman+' > human_reconciledPeaks_intersectAgeSegments.txt')
            pythonInputFiles.append('human_reconciledPeaks_intersectAgeSegments.txt')

            ##########
            # STEP 5:
            # don't need to get reconciled CGIs (like how this step in CGIcentric pipeline gets reconciledPeaks) since these already exist
            # instead use this space to intersect reconciled CGIs with phastCons and age segmentation
             
            # intersect human reconciled peaks with phastCons
            outputString.append('bedtools intersect -wao -a '+ reconciledHuman +' -b '+phastConsHuman+' > human_reconciledCGIs_intersectPhastCons.txt')
            pythonInputFiles.append('human_reconciledCGIs_intersectPhastCons.txt')

            # intersect human reconciled peaks with age segmentation
            outputString.append('bedtools intersect -wao -a '+ reconciledHuman +' -b '+ageSegmentationHuman+' > human_reconciledCGIs_intersectAgeSegments.txt')
            pythonInputFiles.append('human_reconciledCGIs_intersectAgeSegments.txt')


            #STEP 6: run featureCounts to count reads in various peak intervals
            # merge reconciled peaks with full peak lists in each species so that I can calculate RPKMs with only mapped reads 
            outputString.append('bedtools intersect -v -a '+pathToPeaks+speciesA+'_'+tissue+'_'+mark+'_intersection.bed -b speciesA_reconciledPeaks.bed | cat - speciesA_reconciledPeaks.bed > cat_A.bed')
            outputString.append('sort -k1,1 -k2,2n cat_A.bed > sort_A.bed')
            outputString.append('cut -f 1,2,3,4 sort_A.bed > forFeatureCounts_A.bed')

            outputString.append('bedtools intersect -v -a '+pathToPeaks+speciesB+'_'+tissue+'_'+mark+'_intersection.bed -b speciesB_reconciledPeaks.bed | cat - speciesB_reconciledPeaks.bed > cat_B.bed')
            outputString.append('sort -k1,1 -k2,2n cat_B.bed > sort_B.bed')
            outputString.append('cut -f 1,2,3,4 sort_B.bed > forFeatureCounts_B.bed')

            #turn into GTFs
            outputString.append('python /home/ak2267/Scripts/makeGTF_nameFromCol4.py forFeatureCounts_A.bed forFeatureCounts_A.gtf')
            outputString.append('python /home/ak2267/Scripts/makeGTF_nameFromCol4.py forFeatureCounts_B.bed forFeatureCounts_B.gtf')
            outputString.append('module load Subread/2.0.0-GCC-7.3.0-2.30')
            
            #run featureCounts
            A_fileList = []
            B_fileList = []
            for replicate in Rep_Dict[speciesA+'_'+tissue+'_'+mark]:
               outputString.append('featureCounts -t exon -g gene_id -a forFeatureCounts_A.gtf -o speciesA_'+mark+'_'+str(replicate)+'.counts /gpfs/gibbs/pi/noonan/ak2267/Roller/bam/'+speciesA+'_'+tissue+'_'+mark+'_'+str(replicate)+'_bowtie2_filtered.bam')
               A_fileList.append('speciesA_'+mark+'_'+str(replicate)+'.counts')
            for replicate in Rep_Dict[speciesB+'_'+tissue+'_'+mark]:
               outputString.append('featureCounts -t exon -g gene_id -a forFeatureCounts_B.gtf -o speciesB_'+mark+'_'+str(replicate)+'.counts /gpfs/gibbs/pi/noonan/ak2267/Roller/bam/'+speciesB+'_'+tissue+'_'+mark+'_'+str(replicate)+'_bowtie2_filtered.bam')
               B_fileList.append('speciesB_'+mark+'_'+str(replicate)+'.counts')
               
            # generate file with RPM for each reconciledPeak (simpler than doing this in the final python script)
            outputString.append('python /gpfs/gibbs/pi/noonan/ak2267/speciesPairs/Scripts/countRPM.py '+str(len(A_fileList))+' '+' '.join(A_fileList)+' > speciesA_'+mark+'.signal')
            outputString.append('python /gpfs/gibbs/pi/noonan/ak2267/speciesPairs/Scripts/countRPM.py '+str(len(B_fileList))+' '+' '.join(B_fileList)+' > speciesB_'+mark+'.signal')

            pythonInputFiles.append('speciesA_'+mark+'.signal')
            pythonInputFiles.append('speciesB_'+mark+'.signal')

            # STEP 7: run bedtools getfasta and faCount to get CpG numbers in reconciled CGIs, and in reconciled Peaks ***START HERE
            outputString.append('bedtools getfasta -name -fi /gpfs/gibbs/pi/noonan/ak2267/genomes/'+speciesA+'.fa -bed ' +reconciledCGI_A+ ' > speciesA_reconciledCGI.fa ; faCount speciesA_reconciledCGI.fa > speciesA_reconciledCGI.faCount')
            outputString.append('bedtools getfasta -name -fi /gpfs/gibbs/pi/noonan/ak2267/genomes/'+speciesB+'.fa -bed ' +reconciledCGI_B+ ' > speciesB_reconciledCGI.fa ; faCount speciesB_reconciledCGI.fa > speciesB_reconciledCGI.faCount')
            outputString.append('bedtools getfasta -name -fi /gpfs/gibbs/pi/noonan/ak2267/genomes/'+speciesA+'.fa -bed speciesA_reconciledPeaks.bed > speciesA_reconciledPeaks.fa ; faCount speciesA_reconciledPeaks.fa > speciesA_reconciledPeaks.faCount')
            outputString.append('bedtools getfasta -name -fi /gpfs/gibbs/pi/noonan/ak2267/genomes/'+speciesB+'.fa -bed speciesB_reconciledPeaks.bed > speciesB_reconciledPeaks.fa ; faCount speciesB_reconciledPeaks.fa > speciesB_reconciledPeaks.faCount')

            pythonInputFiles.append('speciesA_reconciledCGI.faCount')
            pythonInputFiles.append('speciesB_reconciledCGI.faCount')
            pythonInputFiles.append('speciesA_reconciledPeaks.faCount')
            pythonInputFiles.append('speciesB_reconciledPeaks.faCount')

            # STEP 7.5: append flag for whether to filter based on repeats & lengths or not
            pythonInputFiles.append(filterFlag)

            # STEP 8: if species pairs are possible to consider with phastBias, intersect with phastBias tracts
            #if speciesCombo in phastBiasCombos:
            #    outputString.append('bedtools intersect -wao -a ' + reconciledHuman + ' -b '+pathToPhastBias+phastBias_SpeciesTranslator[speciesA]+'_allChr_phastBias_merged.bed > humanCGI_intersectPhastBias_speciesA.txt')
            #    outputString.append('bedtools intersect -wao -a ' + reconciledHuman + ' -b '+pathToPhastBias+phastBias_SpeciesTranslator[speciesB]+'_allChr_phastBias_merged.bed > humanCGI_intersectPhastBias_speciesB.txt')
            #    outputString.append('bedtools intersect -wao -a human_reconciledPeaks.bed -b '+pathToPhastBias+phastBias_SpeciesTranslator[speciesA]+'_allChr_phastBias_merged.bed > humanPeaks_intersectPhastBias_speciesA.txt')
            #    outputString.append('bedtools intersect -wao -a human_reconciledPeaks.bed -b '+pathToPhastBias+phastBias_SpeciesTranslator[speciesB]+'_allChr_phastBias_merged.bed > humanPeaks_intersectPhastBias_speciesB.txt')

            #    pythonInputFiles.append('humanCGI_intersectPhastBias_speciesA.txt')
            #    pythonInputFiles.append('humanCGI_intersectPhastBias_speciesB.txt')
            #    pythonInputFiles.append('humanPeaks_intersectPhastBias_speciesA.txt')
            #    pythonInputFiles.append('humanPeaks_intersectPhastBias_speciesB.txt')
                
            # STEP 9: python script to analyze the above intersection files & output final counts
            #echo $species $tissue $mark
            outputString.append('python /gpfs/gibbs/pi/noonan/ak2267/speciesPairs/Scripts/countFinalOverlaps_speciesPairs_peakCentric.py '+' '.join(pythonInputFiles)+' >> '+workingDir+'/Roller_summaryFiles/'+speciesCombo+'_'+tissue+'_'+mark+'.txt')

            print('; '.join(outputString))
            
