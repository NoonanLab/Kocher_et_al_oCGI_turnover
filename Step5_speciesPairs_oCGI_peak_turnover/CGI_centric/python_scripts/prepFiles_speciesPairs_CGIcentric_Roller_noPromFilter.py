# 1/10/23
# This script is modified from prepFiles_speciesPairs_CGIcentric_Roller.py
# In the original script, there is a step where oCGIs are filtered out if they overlap a PEAK that intersects a genomic feature
# I want this filtering step in my pipeline for Fig 3-4, but NOT for Fig 2
# This script generates full tables for brain H3K27ac, although the intersection with peaks will not be used downstream for Fig 2

# 7/29/22
# Purpose: make jobFile for doing analysis to make full spreadsheet for speciesPairs, CGIcentric
# Where each line in the spreadsheet is a CGI that maps between the species
# This file only does Roller histone mod data - see other scripts for Liver TF data and Noonan histone mod data
# Modified 9/15/22 to incorporate phastBias

# python prepFiles_speciesPairs_CGIcentric_Roller_noPromFilter.py /gpfs/gibbs/pi/noonan/ak2267/Roller/peaks/repsToUseRoller.txt repeatFilter > 230110_speciesPairs_CGIcentric_Roller_noPromFilter_jobFile.txt

import sys

speciesArray = ['rheMac10','calJac4','mm39','rn7','susScr11','canFam6','felCat9','equCab3']
Tissues = ['brain']
Marks = ['H3K27ac']

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
workingDir = '/gpfs/gibbs/pi/noonan/ak2267/speciesPairs/CGIcentric'
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
            reconciledA = pathToConsensusCGIs+speciesA+'_'+speciesB+'/'+speciesA+'_'+speciesB+'_speciesA_reconciledCGI.bed'
            reconciledB = pathToConsensusCGIs+speciesA+'_'+speciesB+'/'+speciesA+'_'+speciesB+'_speciesB_reconciledCGI.bed'
            reconciledHuman = pathToConsensusCGIs+speciesA+'_'+speciesB+'/'+speciesA+'_'+speciesB+'_hg38_reconciledCGI.bed'
            phastConsHuman = '/home/ak2267/genomes/phastCons/phastConsElements100way_hg38.bed'
            ageSegmentationHuman = '/home/ak2267/project/EnhancerAge/clade_specific_elements_hg38.bed'

            outputString = []

            #########
            # STEP 0: append info on where to work and source files
            outputString.append('cd '+workingDir+' ; source /home/ak2267/.bashrc ; source /home/ak2267/.bash_profile')
            outputString.append('mkdir Roller_'+speciesCombo+'_'+tissue+'_'+mark+' ; cd Roller_'+speciesCombo+'_'+tissue+'_'+mark)

            ######### (steps 1-3 are in reconciledCGI files)
            # STEP 4: intersect reconciled CGIs with several files in each species:
            # peaks for whichever mark to REMOVE CGIs overlapping promoter peaks - include exons, lncRNA, ncRNA, pseudo peaks (so featureAnnotations files 1,2,3,4,5,6)
            # also exclude those overlapping FANTOM peaks and ENCODE blacklist regions as these may not be reliable
            # CGIs to get CGI info
            # rmsk
            # phastCons /home/ak2267/genomes/phastCons
            # age segments
            # peaks for analysis

            pythonInputFiles = []

            # restrict to CGIs that don't overlap a promoter peak in either species
            # ACTUALLY DON'T DO THIS - remove in this version of the script and simply copy files over
            
            outputString.append('cp '+reconciledA+' speciesA_reconciledCGI_noFeaturePeaks.bed')
            outputString.append('cp '+reconciledB+' speciesB_reconciledCGI_noFeaturePeaks.bed')
            outputString.append('cp '+reconciledHuman+' human_reconciledCGI_noFeaturePeaks.bed')
            # for all next steps use speciesA_reconciledCGI_noFeaturePeaks.bed and speciesB_reconciledCGI_noFeaturePeaks.bed

            # intersect each reconciled list with CGI files to get CGI info within python script
            # then: {speciesA}_{speciesB}/{speciesA}_{speciesB}_speciesA_reconciledCGI.bed
            outputString.append('bedtools intersect -wao -a speciesA_reconciledCGI_noFeaturePeaks.bed -b /home/ak2267/genomes/CGI/UCSC_AL/'+speciesA+'_CGIsAL.bed > speciesA_reconciledCGI_intersectSpeciesA_CGI.txt')
            outputString.append('bedtools intersect -wao -a speciesB_reconciledCGI_noFeaturePeaks.bed -b /home/ak2267/genomes/CGI/UCSC_AL/'+speciesB+'_CGIsAL.bed > speciesB_reconciledCGI_intersectSpeciesB_CGI.txt')
            
            pythonInputFiles.append('speciesA_reconciledCGI_intersectSpeciesA_CGI.txt')
            pythonInputFiles.append('speciesB_reconciledCGI_intersectSpeciesB_CGI.txt')

            # intersect each reconciled list with RepeatMasker
            outputString.append('bedtools intersect -wao -a speciesA_reconciledCGI_noFeaturePeaks.bed -b /home/ak2267/genomes/rmsk/rmsk_'+speciesA+'_merged.bed > speciesA_reconciledCGI_intersectSpeciesA_rmsk.txt')
            outputString.append('bedtools intersect -wao -a speciesB_reconciledCGI_noFeaturePeaks.bed -b /home/ak2267/genomes/rmsk/rmsk_'+speciesB+'_merged.bed > speciesB_reconciledCGI_intersectSpeciesB_rmsk.txt')

            pythonInputFiles.append('speciesA_reconciledCGI_intersectSpeciesA_rmsk.txt')
            pythonInputFiles.append('speciesB_reconciledCGI_intersectSpeciesB_rmsk.txt')

            # intersect human reconciled list with phastCons elements
            outputString.append('bedtools intersect -wao -a human_reconciledCGI_noFeaturePeaks.bed -b '+phastConsHuman+' > human_reconciledCGI_intersectPhastCons.txt')
            pythonInputFiles.append('human_reconciledCGI_intersectPhastCons.txt')

            # intersect human reconcile list with age segmentation map
            outputString.append('bedtools intersect -wao -a human_reconciledCGI_noFeaturePeaks.bed -b '+ageSegmentationHuman+' > human_reconciledCGI_intersectAgeSegments.txt')
            pythonInputFiles.append('human_reconciledCGI_intersectAgeSegments.txt')

            ##########
            # STEP 5:
            # run same pipeline as for getting reconciled CGIs to get orthologous peak regions in each species

            # intersect each reconciled list with mark in the appropriate species & write complete output
            # so that mark overlap is captured even if reconciledPeaks fails
            outputString.append('bedtools intersect -wao -a speciesA_reconciledCGI_noFeaturePeaks.bed -b '+pathToPeaks+speciesA+'_'+tissue+'_'+mark+'_intersection.bed > speciesA_reconciledCGI_intersect'+mark+'.txt')
            outputString.append('bedtools intersect -wao -a speciesB_reconciledCGI_noFeaturePeaks.bed -b '+pathToPeaks+speciesB+'_'+tissue+'_'+mark+'_intersection.bed > speciesB_reconciledCGI_intersect'+mark+'.txt')
            pythonInputFiles.append('speciesA_reconciledCGI_intersect'+mark+'.txt')
            pythonInputFiles.append('speciesB_reconciledCGI_intersect'+mark+'.txt')

            # begin pipeline for reconciledPeaks
            # restrict to peaks overlapping reconciled CGI and put name in column 4
            outputString.append('bedtools intersect -wa -u -a '+pathToPeaks+speciesA+'_'+tissue+'_'+mark+'_intersection.bed -b speciesA_reconciledCGI_noFeaturePeaks.bed > speciesA_peaksWithCGI.bed')
            outputString.append('bedtools intersect -wa -u -a '+pathToPeaks+speciesB+'_'+tissue+'_'+mark+'_intersection.bed -b speciesB_reconciledCGI_noFeaturePeaks.bed > speciesB_peaksWithCGI.bed')
            outputString.append('awk \'{ print $1\"\\t\"$2\"\\t\"$3\"\\ta_\"NR }\' speciesA_peaksWithCGI.bed > speciesA_peaksWithCGI_4col.bed')
            outputString.append('awk \'{ print $1\"\\t\"$2\"\\t\"$3\"\\tb_\"NR }\' speciesB_peaksWithCGI.bed > speciesB_peaksWithCGI_4col.bed')

            # lift to human
            chainAtoHg38 = '/home/ak2267/genomes/chain/'+speciesA+'ToHg38.over.chain.gz'
            chainBtoHg38 = '/home/ak2267/genomes/chain/'+speciesB+'ToHg38.over.chain.gz'
            chainHg38toA = '/home/ak2267/genomes/chain/hg38To'+speciesA[0].upper()+speciesA[1:]+'.over.chain.gz'
            chainHg38toB = '/home/ak2267/genomes/chain/hg38To'+speciesB[0].upper()+speciesB[1:]+'.over.chain.gz'

            outputString.append('liftOver -minMatch=0.3 speciesA_peaksWithCGI_4col.bed '+chainAtoHg38+' speciesA_peaksWithCGI_LOhuman.bed unMapped')
            outputString.append('liftOver -minMatch=0.3 speciesB_peaksWithCGI_4col.bed '+chainBtoHg38+' speciesB_peaksWithCGI_LOhuman.bed unMapped')

            # lift back to other species and keep sites in human that lift to the same peak as the one they started as (analagous to piepline in consensusCGIs)
            outputString.append('liftOver -minMatch=0.3 speciesA_peaksWithCGI_LOhuman.bed '+chainHg38toA+' speciesA_peaksWithCGI_liftedBackToA.bed unMapped')
            outputString.append('liftOver -minMatch=0.3 speciesB_peaksWithCGI_LOhuman.bed '+chainHg38toB+' speciesB_peaksWithCGI_liftedBackToB.bed unMapped')
            outputString.append('bedtools intersect -wao -a speciesA_peaksWithCGI_liftedBackToA.bed -b speciesA_peaksWithCGI_4col.bed > speciesA_checkMappingBack.txt')
            outputString.append('bedtools intersect -wao -a speciesB_peaksWithCGI_liftedBackToB.bed -b speciesB_peaksWithCGI_4col.bed > speciesB_checkMappingBack.txt')
            outputString.append('python /gpfs/gibbs/pi/noonan/ak2267/speciesPairs/consensusCGIs/Scripts/restrictToSitesThatMapBack.py speciesA_checkMappingBack.txt > speciesA_peaksWithCGI_mapBack_inAcoord.bed')
            outputString.append('python /gpfs/gibbs/pi/noonan/ak2267/speciesPairs/consensusCGIs/Scripts/restrictToSitesThatMapBack.py speciesB_checkMappingBack.txt > speciesB_peaksWithCGI_mapBack_inBcoord.bed')
            outputString.append('python /gpfs/gibbs/pi/noonan/ak2267/speciesPairs/consensusCGIs/Scripts/restrictToLO.py speciesA_peaksWithCGI_LOhuman.bed speciesA_peaksWithCGI_mapBack_inAcoord.bed > speciesA_peaksWithCGI_LOhuman_mapBack.bed')
            outputString.append('python /gpfs/gibbs/pi/noonan/ak2267/speciesPairs/consensusCGIs/Scripts/restrictToLO.py speciesB_peaksWithCGI_LOhuman.bed speciesB_peaksWithCGI_mapBack_inBcoord.bed > speciesB_peaksWithCGI_LOhuman_mapBack.bed')
               
            # merge, remove features in human, sort file names to be a_# then b_#
            outputString.append('cat speciesA_peaksWithCGI_LOhuman_mapBack.bed speciesB_peaksWithCGI_LOhuman_mapBack.bed > cat.bed')
            outputString.append('sort -k1,1 -k2,2n cat.bed > sort.bed')
            outputString.append('bedtools merge -c 4 -o collapse -i sort.bed > mergePeaks_human.bed')
            outputString.append('python /gpfs/gibbs/pi/noonan/ak2267/speciesPairs/consensusCGIs/Scripts/sortRegionNames.py mergePeaks_human.bed > mergePeaks_human_sortedNames.bed')

            # lift back out to speciesA and speciesB
            outputString.append('liftOver -minMatch=0.3 mergePeaks_human_sortedNames.bed '+chainHg38toA+' mergePeaks_liftToA.bed unMapped')
            outputString.append('liftOver -minMatch=0.3 mergePeaks_human_sortedNames.bed '+chainHg38toB+' mergePeaks_liftToB.bed unMapped')

            # lift back to human and restrict list in speciesA and speciesB to only those that lift back to the same CGI in human (analagous to above)
            outputString.append('liftOver -minMatch=0.3 mergePeaks_liftToA.bed '+chainAtoHg38+' mergePeaks_liftToA_backToHuman.bed unMapped')
            outputString.append('liftOver -minMatch=0.3 mergePeaks_liftToB.bed '+chainBtoHg38+' mergePeaks_liftToB_backToHuman.bed unMapped')
            outputString.append('bedtools intersect -wao -a mergePeaks_liftToA_backToHuman.bed -b mergePeaks_human_sortedNames.bed > speciesA_checkMappingBackInHuman.txt')
            outputString.append('bedtools intersect -wao -a mergePeaks_liftToB_backToHuman.bed -b mergePeaks_human_sortedNames.bed > speciesB_checkMappingBackInHuman.txt')
            outputString.append('python /gpfs/gibbs/pi/noonan/ak2267/speciesPairs/consensusCGIs/Scripts/restrictToSitesThatMapBack.py speciesA_checkMappingBackInHuman.txt > speciesA_mergePeaks_thatMapBackToHuman_inHuman.bed')
            outputString.append('python /gpfs/gibbs/pi/noonan/ak2267/speciesPairs/consensusCGIs/Scripts/restrictToSitesThatMapBack.py speciesB_checkMappingBackInHuman.txt > speciesB_mergePeaks_thatMapBackToHuman_inHuman.bed')
            outputString.append('python /gpfs/gibbs/pi/noonan/ak2267/speciesPairs/consensusCGIs/Scripts/restrictToLO.py mergePeaks_liftToA.bed speciesA_mergePeaks_thatMapBackToHuman_inHuman.bed > mergePeaks_speciesA_mapBack.bed')
            outputString.append('python /gpfs/gibbs/pi/noonan/ak2267/speciesPairs/consensusCGIs/Scripts/restrictToLO.py mergePeaks_liftToB.bed speciesB_mergePeaks_thatMapBackToHuman_inHuman.bed > mergePeaks_speciesB_mapBack.bed')

            # restrict lists (in A, in B) to those present in both A and B - also do for human list
            outputString.append('python /gpfs/gibbs/pi/noonan/ak2267/speciesPairs/consensusCGIs/Scripts/restrictToLO.py mergePeaks_speciesA_mapBack.bed mergePeaks_speciesB_mapBack.bed > speciesA_reconciledPeaks.bed')
            outputString.append('python /gpfs/gibbs/pi/noonan/ak2267/speciesPairs/consensusCGIs/Scripts/restrictToLO.py mergePeaks_speciesB_mapBack.bed mergePeaks_speciesA_mapBack.bed > speciesB_reconciledPeaks.bed')
            outputString.append('python /gpfs/gibbs/pi/noonan/ak2267/speciesPairs/consensusCGIs/Scripts/restrictToLO.py mergePeaks_human_sortedNames.bed mergePeaks_speciesA_mapBack.bed > human_intermediate_reconciledPeaks.bed')
            outputString.append('python /gpfs/gibbs/pi/noonan/ak2267/speciesPairs/consensusCGIs/Scripts/restrictToLO.py human_intermediate_reconciledPeaks.bed mergePeaks_speciesB_mapBack.bed > human_reconciledPeaks.bed')
            
            # intersect human reconciled peaks with phastCons
            outputString.append('bedtools intersect -wao -a human_reconciledPeaks.bed -b '+phastConsHuman+' > human_reconciledPeaks_intersectPhastCons.txt')
            pythonInputFiles.append('human_reconciledPeaks_intersectPhastCons.txt')

            # intersect human reconciled peaks with age segmentation
            outputString.append('bedtools intersect -wao -a human_reconciledPeaks.bed -b '+ageSegmentationHuman+' > human_reconciledPeaks_intersectAgeSegments.txt')
            pythonInputFiles.append('human_reconciledPeaks_intersectAgeSegments.txt')

            # intersect reconciled CGIs with reconciled Peaks
            outputString.append('bedtools intersect -wao -a speciesA_reconciledCGI_noFeaturePeaks.bed -b speciesA_reconciledPeaks.bed > speciesA_reconciledCGI_intersect_reconciledPeaks.txt')
            outputString.append('bedtools intersect -wao -a speciesB_reconciledCGI_noFeaturePeaks.bed -b speciesB_reconciledPeaks.bed > speciesB_reconciledCGI_intersect_reconciledPeaks.txt')
            pythonInputFiles.append('speciesA_reconciledCGI_intersect_reconciledPeaks.txt')
            pythonInputFiles.append('speciesB_reconciledCGI_intersect_reconciledPeaks.txt')

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

            # STEP 7: run bedtools getfasta and faCount to get CpG numbers in reconciled CGIs, and in reconciled Peaks
            outputString.append('bedtools getfasta -name -fi /gpfs/gibbs/pi/noonan/ak2267/genomes/'+speciesA+'.fa -bed speciesA_reconciledCGI_noFeaturePeaks.bed > speciesA_reconciledCGI.fa ; faCount speciesA_reconciledCGI.fa > speciesA_reconciledCGI.faCount')
            outputString.append('bedtools getfasta -name -fi /gpfs/gibbs/pi/noonan/ak2267/genomes/'+speciesB+'.fa -bed speciesB_reconciledCGI_noFeaturePeaks.bed > speciesB_reconciledCGI.fa ; faCount speciesB_reconciledCGI.fa > speciesB_reconciledCGI.faCount')
            outputString.append('bedtools getfasta -name -fi /gpfs/gibbs/pi/noonan/ak2267/genomes/'+speciesA+'.fa -bed speciesA_reconciledPeaks.bed > speciesA_reconciledPeaks.fa ; faCount speciesA_reconciledPeaks.fa > speciesA_reconciledPeaks.faCount')
            outputString.append('bedtools getfasta -name -fi /gpfs/gibbs/pi/noonan/ak2267/genomes/'+speciesB+'.fa -bed speciesB_reconciledPeaks.bed > speciesB_reconciledPeaks.fa ; faCount speciesB_reconciledPeaks.fa > speciesB_reconciledPeaks.faCount')

            pythonInputFiles.append('speciesA_reconciledCGI.faCount')
            pythonInputFiles.append('speciesB_reconciledCGI.faCount')
            pythonInputFiles.append('speciesA_reconciledPeaks.faCount')
            pythonInputFiles.append('speciesB_reconciledPeaks.faCount')

            # STEP 7.5: append flag for whether to filter based on repeats & lengths or not
            pythonInputFiles.append(filterFlag)

            # STEP 8: if species pairs are possible to consider with phastBias, intersect with phastBias tracts
            if speciesCombo in phastBiasCombos:
                outputString.append('bedtools intersect -wao -a human_reconciledCGI_noFeaturePeaks.bed -b '+pathToPhastBias+phastBias_SpeciesTranslator[speciesA]+'_allChr_phastBias_merged.bed > humanCGI_intersectPhastBias_speciesA.txt')
                outputString.append('bedtools intersect -wao -a human_reconciledCGI_noFeaturePeaks.bed -b '+pathToPhastBias+phastBias_SpeciesTranslator[speciesB]+'_allChr_phastBias_merged.bed > humanCGI_intersectPhastBias_speciesB.txt')
                outputString.append('bedtools intersect -wao -a human_reconciledPeaks.bed -b '+pathToPhastBias+phastBias_SpeciesTranslator[speciesA]+'_allChr_phastBias_merged.bed > humanPeaks_intersectPhastBias_speciesA.txt')
                outputString.append('bedtools intersect -wao -a human_reconciledPeaks.bed -b '+pathToPhastBias+phastBias_SpeciesTranslator[speciesB]+'_allChr_phastBias_merged.bed > humanPeaks_intersectPhastBias_speciesB.txt')

                pythonInputFiles.append('humanCGI_intersectPhastBias_speciesA.txt')
                pythonInputFiles.append('humanCGI_intersectPhastBias_speciesB.txt')
                pythonInputFiles.append('humanPeaks_intersectPhastBias_speciesA.txt')
                pythonInputFiles.append('humanPeaks_intersectPhastBias_speciesB.txt')
                
            # STEP 9: python script to analyze the above intersection files & output final counts
            #echo $species $tissue $mark
            outputString.append('python /gpfs/gibbs/pi/noonan/ak2267/speciesPairs/Scripts/countFinalOverlaps_speciesPairs_CGIcentric.py '+' '.join(pythonInputFiles)+' >> '+workingDir+'/SummaryFiles_noPromFilter/'+speciesCombo+'_'+tissue+'_'+mark+'.txt')

            print('; '.join(outputString))
            
