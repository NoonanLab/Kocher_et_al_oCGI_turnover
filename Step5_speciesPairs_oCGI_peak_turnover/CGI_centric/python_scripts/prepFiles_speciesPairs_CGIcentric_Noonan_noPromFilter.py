# 1/10/23
# This script is modified from prepFiles_speciesPairs_CGIcentric_Roller.py
# In the original script, there is a step where oCGIs are filtered out if they overlap a PEAK that intersects a genomic feature
# I want this filtering step in my pipeline for Fig 3-4, but NOT for Fig 2
# This script generates full tables for brain H3K27ac at timePoint 1, although the intersection with peaks will not be used downstream for Fig 2

# 8/10/22
# Purpose: make jobFile for doing analysis to make full spreadsheet for speciesPairs, CGIcentric
# Where each line in the spreadsheet is a CGI that maps between the species
# This file only Noonan data (brain and limb) - see other scripts for Roller histone mod data and Liver TF data
# Modified 9/15/22 to incorporate phastBias
# Modified 9/26/22 to fix issue with some samples not getting bigWig data incorporated - bad file paths

import sys

Tissues = ['brain']
speciesArray = ['hg19','rheMac2','mm9']

Marks_Dict = {}
Marks_Dict['brain'] = ['ac']

phastBiasCombos = ['hg19_rheMac2']
phastBias_SpeciesTranslator = {'hg19':'hg38inhg19', 'rheMac2':'rheMac8inhg19'}

filterFlag = str(sys.argv[1]) # repeatFilter vs noRepeatFilter

# store path to phastBias
pathToPhastBias = '/gpfs/gibbs/pi/noonan/ak2267/speciesPairs/SequenceAnalysis/gBGC_analysis/phastBias_bed/'
# then: ${species}_allChr_phastBias_merged.bed - note this is in hg38 coordinates

# store timePoints for comparison in dictionary:
# PEAK TIMEPOINTS AND BIG WIG TIMEPOINTS DO NOT ALWAYS MATCH - see mouse 17F/17O vs e17F/e17O
TimePoints_Dict = {}
for tissue in ['brain']:
    TimePoints_Dict[tissue] = {}
TimePoints_Dict['brain']['hg19'] = ['CS16','CS23','F2F','F2O']
TimePoints_Dict['brain']['rheMac2'] = ['.','e55','e79F','e79O']
TimePoints_Dict['brain']['mm9'] =['e11','e14','17F','17O']

# store paths to peak files in dictionary - all require timePoint to be inserted between 0 and 1 in array
PeakPath_Dict = {}
for tissue in ['brain','limb']:
    PeakPath_Dict[tissue] = {}
    for species in ['hg19','rheMac2','mm9']:
        PeakPath_Dict[tissue][species] = {}
PeakPath_Dict['brain']['hg19']['ac'] = ['/home/ak2267/project/EnhancerClasses/hg19/ac/merge_','_overlap_named.bed']
PeakPath_Dict['brain']['rheMac2']['ac'] = ['/home/ak2267/project/EnhancerClasses/rheMac2/ac/merge_','_overlap_rh_named.bed']
PeakPath_Dict['brain']['mm9']['ac'] = ['/home/ak2267/project/EnhancerClasses/mm9/ac/merge_','_overlap_mm_named.bed']

# store paths to bigWig files in dictionary - all require timePoint to be inserted between 0 and 1 in array, and rep between 1 and 2
BigWigPath_Dict = {}
for tissue in ['brain']:
    BigWigPath_Dict[tissue] = {}
    for species in ['hg19','rheMac2','mm9']:
        BigWigPath_Dict[tissue][species] = {}
BigWigPath_Dict['brain']['hg19']['ac'] = ['/gpfs/gibbs/pi/noonan/ak2267/NoonanBrain/hg19/ac/','_ac_rep','.bw']
BigWigPath_Dict['brain']['rheMac2']['ac'] = ['/gpfs/gibbs/pi/noonan/ak2267/NoonanBrain/rheMac2/ac/','_ac_rep','.bw']
BigWigPath_Dict['brain']['mm9']['ac'] = ['/gpfs/gibbs/pi/noonan/ak2267/NoonanBrain/mm9/ac/','_ac_rep','.bw']

# store info on rep number for each combo
RepNumber_Dict = {}
for tissue in ['brain']:
    RepNumber_Dict[tissue] = {}
RepNumber_Dict['brain']['hg19'] = [2,2,2,2]
RepNumber_Dict['brain']['rheMac2'] = [0,1,2,2]
RepNumber_Dict['brain']['mm9'] =[2,2,2,2]

## START PIPELINE

# store filePaths and fileEnds
workingDir = '/gpfs/gibbs/pi/noonan/ak2267/speciesPairs/CGIcentric'
pathToConsensusCGIs = '/gpfs/gibbs/pi/noonan/ak2267/speciesPairs/consensusCGIs/'
# then: {speciesA}_{speciesB}/{speciesA}_{speciesB}_speciesA_reconciledCGI.bed
#    or {speciesA}_{speciesB}/{speciesA}_{speciesB}_speciesB_reconciledCGI.bed
pathToAnnotations = '/home/ak2267/genomes/RefSeq/featureAnnotations/'
# then: ${species}_allFeatures.bed
pathToFANTOM = '/home/ak2267/genomes/FANTOM/'
# then: FANTOM_TSS_${species}.bed
pathToBlacklist = '/home/ak2267/genomes/blacklist/'
# then: blacklist_${species}.bed

SpeciesCombos = ['hg19_rheMac2','hg19_mm9','rheMac2_mm9']

for speciesCombo in SpeciesCombos:
    for tissue in Tissues:
        for mark in Marks_Dict[tissue]:
            for timePointIndex in [1]:
                #print(speciesCombo+'\t'+tissue+'\t'+mark)
                speciesA = speciesCombo.split('_')[0]
                speciesB = speciesCombo.split('_')[1]

                timePointA = TimePoints_Dict[tissue][speciesA][timePointIndex]
                timePointB = TimePoints_Dict[tissue][speciesB][timePointIndex]

                if timePointA != '.' and timePointB != '.':
                
                    reconciledA = pathToConsensusCGIs+speciesA+'_'+speciesB+'/'+speciesA+'_'+speciesB+'_speciesA_reconciledCGI.bed'
                    reconciledB = pathToConsensusCGIs+speciesA+'_'+speciesB+'/'+speciesA+'_'+speciesB+'_speciesB_reconciledCGI.bed'
                    reconciledHuman = pathToConsensusCGIs+speciesA+'_'+speciesB+'/'+speciesA+'_'+speciesB+'_hg19_reconciledCGI.bed'
                    phastConsHuman = '/home/ak2267/genomes/phastCons/phastConsElements100way_hg19.bed'
                    ageSegmentationHuman = '/home/ak2267/project/EnhancerAge/clade_specific_elements_hg19.bed'

                    pathToPeaksA = timePointA.join(PeakPath_Dict[tissue][speciesA][mark])
                    pathToPeaksB = timePointB.join(PeakPath_Dict[tissue][speciesB][mark])
                    # these are still arrays because they need both timePoint and rep filled in
                    pathToBigWigA = BigWigPath_Dict[tissue][speciesA][mark]
                    pathToBigWigB = BigWigPath_Dict[tissue][speciesB][mark]

                    outputString = []

                    #########
                    # STEP 0: append info on where to work and source files
                    outputString.append('cd '+workingDir+' ; source /home/ak2267/.bashrc ; source /home/ak2267/.bash_profile')
                    outputString.append('mkdir Noonan_'+speciesCombo+'_'+tissue+'_'+mark+'_'+str(timePointIndex)+' ; cd Noonan_'+speciesCombo+'_'+tissue+'_'+mark+'_'+str(timePointIndex))

                    ######### (steps 1-3 are in reconciledCGI files)
                    # STEP 4: intersect reconciled CGIs with several files in each species:
                    # peaks for whichever mark to REMOVE CGIs overlapping promoter peaks - include exons, lncRNA, ncRNA, pseudo peaks (so featureAnnotations files 1,2,3,4,5,6)
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
                    # have to consider that hg19 can be speciesA and doesn't need to be lifted to hg19

                    # intersect each reconciled list with mark in the appropriate species & write complete output
                    # so that mark overlap is captured even if reconciledPeaks fails
                    outputString.append('bedtools intersect -wao -a speciesA_reconciledCGI_noFeaturePeaks.bed -b '+pathToPeaksA+' > speciesA_reconciledCGI_intersect'+mark+'.txt')
                    outputString.append('bedtools intersect -wao -a speciesB_reconciledCGI_noFeaturePeaks.bed -b '+pathToPeaksB+' > speciesB_reconciledCGI_intersect'+mark+'.txt')
                    pythonInputFiles.append('speciesA_reconciledCGI_intersect'+mark+'.txt')
                    pythonInputFiles.append('speciesB_reconciledCGI_intersect'+mark+'.txt')

                    # begin pipeline for reconciledPeaks
                    # restrict to peaks overlapping reconciled CGI and put name in column 4
                    outputString.append('bedtools intersect -wa -u -a '+pathToPeaksA+' -b speciesA_reconciledCGI_noFeaturePeaks.bed > speciesA_peaksWithCGI.bed')
                    outputString.append('bedtools intersect -wa -u -a '+pathToPeaksB+' -b speciesB_reconciledCGI_noFeaturePeaks.bed > speciesB_peaksWithCGI.bed')
                    outputString.append('awk \'{ print $1\"\\t\"$2\"\\t\"$3\"\\ta_\"NR }\' speciesA_peaksWithCGI.bed > speciesA_peaksWithCGI_4col.bed')
                    outputString.append('awk \'{ print $1\"\\t\"$2\"\\t\"$3\"\\tb_\"NR }\' speciesB_peaksWithCGI.bed > speciesB_peaksWithCGI_4col.bed')

                    # lift to human
                    chainAtoHg19 = '/home/ak2267/genomes/chain/'+speciesA+'ToHg19.over.chain.gz'
                    chainBtoHg19 = '/home/ak2267/genomes/chain/'+speciesB+'ToHg19.over.chain.gz'
                    chainHg19toA = '/home/ak2267/genomes/chain/hg19To'+speciesA[0].upper()+speciesA[1:]+'.over.chain.gz'
                    chainHg19toB = '/home/ak2267/genomes/chain/hg19To'+speciesB[0].upper()+speciesB[1:]+'.over.chain.gz'

                    if speciesA == 'hg19':
                        outputString.append('cp speciesA_peaksWithCGI_4col.bed speciesA_peaksWithCGI_LOhuman.bed')
                        outputString.append('liftOver -minMatch=0.3 speciesB_peaksWithCGI_4col.bed '+chainBtoHg19+' speciesB_peaksWithCGI_LOhuman.bed unMapped')
                    if speciesA != 'hg19':
                        outputString.append('liftOver -minMatch=0.3 speciesA_peaksWithCGI_4col.bed '+chainAtoHg19+' speciesA_peaksWithCGI_LOhuman.bed unMapped')
                        outputString.append('liftOver -minMatch=0.3 speciesB_peaksWithCGI_4col.bed '+chainBtoHg19+' speciesB_peaksWithCGI_LOhuman.bed unMapped')

                    # lift back to other species and keep sites in human that lift to the same peak as the one they started as (analagous to piepline in consensusCGIs)
                    if speciesA == 'hg19':
                        outputString.append('cp speciesA_peaksWithCGI_LOhuman.bed speciesA_peaksWithCGI_liftedBackToA.bed')
                        outputString.append('liftOver -minMatch=0.3 speciesB_peaksWithCGI_LOhuman.bed '+chainHg19toB+' speciesB_peaksWithCGI_liftedBackToB.bed unMapped')
                    if speciesA != 'hg19':
                        outputString.append('liftOver -minMatch=0.3 speciesA_peaksWithCGI_LOhuman.bed '+chainHg19toA+' speciesA_peaksWithCGI_liftedBackToA.bed unMapped')
                        outputString.append('liftOver -minMatch=0.3 speciesB_peaksWithCGI_LOhuman.bed '+chainHg19toB+' speciesB_peaksWithCGI_liftedBackToB.bed unMapped')
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
                    if speciesA == 'hg19':
                        outputString.append('cp mergePeaks_human_sortedNames.bed mergePeaks_liftToA.bed')
                        outputString.append('liftOver -minMatch=0.3 mergePeaks_human_sortedNames.bed '+chainHg19toB+' mergePeaks_liftToB.bed unMapped')
                    if speciesA != 'hg19':
                        outputString.append('liftOver -minMatch=0.3 mergePeaks_human_sortedNames.bed '+chainHg19toA+' mergePeaks_liftToA.bed unMapped')
                        outputString.append('liftOver -minMatch=0.3 mergePeaks_human_sortedNames.bed '+chainHg19toB+' mergePeaks_liftToB.bed unMapped')

                    # lift back to human and restrict list in speciesA and speciesB to only those that lift back to the same CGI in human (analagous to above)
                    if speciesA == 'hg19':
                        outputString.append('cp mergePeaks_liftToA.bed mergePeaks_liftToA_backToHuman.bed')
                        outputString.append('liftOver -minMatch=0.3 mergePeaks_liftToB.bed '+chainBtoHg19+' mergePeaks_liftToB_backToHuman.bed unMapped')
                    if speciesA != 'hg19':
                        outputString.append('liftOver -minMatch=0.3 mergePeaks_liftToA.bed '+chainAtoHg19+' mergePeaks_liftToA_backToHuman.bed unMapped')
                        outputString.append('liftOver -minMatch=0.3 mergePeaks_liftToB.bed '+chainBtoHg19+' mergePeaks_liftToB_backToHuman.bed unMapped')
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

                    #STEP 6: run bigWigAverageOverBed to get average signal in peak intervals - THIS DIFFERS FROM ROLLER PIPELINE
                    # merge reconciled peaks with full peak lists in each species so that I can calculate RPKMs with only mapped reads 
                    outputString.append('bedtools intersect -v -a '+pathToPeaksA+' -b speciesA_reconciledPeaks.bed | cat - speciesA_reconciledPeaks.bed > cat_A.bed')
                    outputString.append('sort -k1,1 -k2,2n cat_A.bed > sort_A.bed')
                    outputString.append('cut -f 1,2,3,4 sort_A.bed > forAvgOverBed_A.bed')

                    outputString.append('bedtools intersect -v -a '+pathToPeaksB+' -b speciesB_reconciledPeaks.bed | cat - speciesB_reconciledPeaks.bed > cat_B.bed')
                    outputString.append('sort -k1,1 -k2,2n cat_B.bed > sort_B.bed')
                    outputString.append('cut -f 1,2,3,4 sort_B.bed > forAvgOverBed_B.bed')

                    #run bigWigAverageOverBed
                    A_fileList = []
                    B_fileList = []

                    #pathToPeaksA = timePointA.join(PeakPath_Dict[tissue][speciesA][mark])
                    #pathToPeaksB = timePointB.join(PeakPath_Dict[tissue][speciesB][mark])
                    # these are still arrays because they need both timePoint and rep filled in
                    #pathToBigWigA = BigWigPath_Dict[tissue][speciesA][mark]
                    #pathToBigWigB = BigWigPath_Dict[tissue][speciesB][mark]

                    #PeakPath_Dict['limb']['mm9']['ac'] = ['/home/ak2267/project/EnhancerClasses/Limb/mm9/merge_','_overlap_named.bed']
                    #BigWigPath_Dict['brain']['hg19']['ac'] = ['/gpfs/gibbs/pi/noonan/ak2267/NoonanBrain/hg19/ac/','_ac_rep','.bw']
                    
                    # change mouse timePoint name (bigWigs are named differently than peak files)
                    fixMouseTimePoints = {'17F':'e17F','17O':'e17O'}
                    if timePointB in fixMouseTimePoints:
                    	timePointB = fixMouseTimePoints[timePointB]

					# add command for speciesA 
					# naming is different for this replicate of rhesus - fix manually
                    if timePointIndex == 1 and tissue == 'limb' and speciesA == 'rheMac2':
                        outputString.append('bigWigAverageOverBed '+pathToBigWigA[0]+'e31'+pathToBigWigA[1]+'1'+pathToBigWigA[2]+' forAvgOverBed_A.bed speciesA_1.counts')
                        outputString.append('bigWigAverageOverBed '+pathToBigWigA[0]+'e32'+pathToBigWigA[1]+'1'+pathToBigWigA[2]+' forAvgOverBed_A.bed speciesA_2.counts')
                        outputString.append('bigWigAverageOverBed '+pathToBigWigA[0]+'e33'+pathToBigWigA[1]+'1'+pathToBigWigA[2]+' forAvgOverBed_A.bed speciesA_3.counts')
                        A_fileList.append('speciesA_1.counts')
                        A_fileList.append('speciesA_2.counts')
                        A_fileList.append('speciesA_3.counts')

                    else:
                        for replicate in range(1,RepNumber_Dict[tissue][speciesA][timePointIndex]+1):
                            outputString.append('bigWigAverageOverBed '+pathToBigWigA[0]+timePointA+pathToBigWigA[1]+str(replicate)+pathToBigWigA[2]+' forAvgOverBed_A.bed speciesA_'+str(replicate)+'.counts')
                            A_fileList.append('speciesA_'+str(replicate)+'.counts')
                        
                    # add command for speciesB   
                    if timePointIndex == 1 and tissue == 'limb' and speciesB == 'rheMac2':
                        outputString.append('bigWigAverageOverBed '+pathToBigWigB[0]+'e31'+pathToBigWigB[1]+'1'+pathToBigWigB[2]+' forAvgOverBed_B.bed speciesB_1.counts')
                        outputString.append('bigWigAverageOverBed '+pathToBigWigB[0]+'e32'+pathToBigWigB[1]+'1'+pathToBigWigB[2]+' forAvgOverBed_B.bed speciesB_2.counts')
                        outputString.append('bigWigAverageOverBed '+pathToBigWigB[0]+'e33'+pathToBigWigB[1]+'1'+pathToBigWigB[2]+' forAvgOverBed_B.bed speciesB_3.counts')
                        B_fileList.append('speciesB_1.counts')
                        B_fileList.append('speciesB_2.counts')
                        B_fileList.append('speciesB_3.counts')

                    else:
                        for replicate in range(1,RepNumber_Dict[tissue][speciesB][timePointIndex]+1):
                            outputString.append('bigWigAverageOverBed '+pathToBigWigB[0]+timePointB+pathToBigWigB[1]+str(replicate)+pathToBigWigB[2]+' forAvgOverBed_B.bed speciesB_'+str(replicate)+'.counts')
                            B_fileList.append('speciesB_'+str(replicate)+'.counts')

                    # generate file with averageSignal for each reconciledPeak (simpler than doing this in the final python script)
                    outputString.append('python /gpfs/gibbs/pi/noonan/ak2267/speciesPairs/Scripts/countBigWigSignal.py '+str(len(A_fileList))+' '+' '.join(A_fileList)+' > speciesA_'+mark+'.signal')
                    outputString.append('python /gpfs/gibbs/pi/noonan/ak2267/speciesPairs/Scripts/countBigWigSignal.py '+str(len(B_fileList))+' '+' '.join(B_fileList)+' > speciesB_'+mark+'.signal')

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
                    # USING ROLLER SCRIPT BUT NOTE THAT THE COLUMNS TITLED RPM AND RPKM NEED TO BE INTERPRETED DIFFERENTLY
                    outputString.append('python /gpfs/gibbs/pi/noonan/ak2267/speciesPairs/Scripts/countFinalOverlaps_speciesPairs_CGIcentric.py '+' '.join(pythonInputFiles)+' >> '+workingDir+'/SummaryFiles_noPromFilter/'+speciesCombo+'_'+tissue+'_'+mark+'_'+str(timePointIndex)+'.txt')

                    print('; '.join(outputString))

