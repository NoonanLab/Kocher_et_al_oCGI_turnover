# 11/30/22
# Purpose: make jobFile for doing analysis to make full spreadsheet for speciesPairs, peakCentric
# Where each line in the spreadsheet is a PEAK that maps between the species
# Modified 12/2/22 to process Noonan data - see other file for Roller data (prepFiles_speciesPairs_peakCentric_Roller.py)

import sys

Tissues = ['brain','limb']
speciesArray = ['hg19','rheMac2','mm9']

Marks_Dict = {}
Marks_Dict['brain'] = ['ac','me2']
Marks_Dict['limb'] = ['ac']

filterFlag = str(sys.argv[1]) # repeatFilter vs noRepeatFilter

# store timePoints for comparison in dictionary:
# PEAK TIMEPOINTS AND BIG WIG TIMEPOINTS DO NOT ALWAYS MATCH - see mouse 17F/17O vs e17F/e17O
TimePoints_Dict = {}
for tissue in ['brain','limb']:
    TimePoints_Dict[tissue] = {}
TimePoints_Dict['brain']['hg19'] = ['CS16','CS23','F2F','F2O']
TimePoints_Dict['brain']['rheMac2'] = ['.','e55','e79F','e79O']
TimePoints_Dict['brain']['mm9'] =['e11','e14','17F','17O']
TimePoints_Dict['limb']['hg19'] = ['E33','E41','E44','E47']
TimePoints_Dict['limb']['rheMac2'] = ['.','e31','e36','.']
TimePoints_Dict['limb']['mm9'] =['e10.5','e11.5','e12.5','e13.5']

# store paths to peak files in dictionary - all require timePoint to be inserted between 0 and 1 in array
PeakPath_Dict = {}
for tissue in ['brain','limb']:
    PeakPath_Dict[tissue] = {}
    for species in ['hg19','rheMac2','mm9']:
        PeakPath_Dict[tissue][species] = {}
PeakPath_Dict['brain']['hg19']['ac'] = ['/home/ak2267/project/EnhancerClasses/hg19/ac/merge_','_overlap_named.bed']
PeakPath_Dict['brain']['hg19']['me2'] = ['/home/ak2267/project/EnhancerClasses/hg19/me2/merge_','_me2_overlap_named.bed']
PeakPath_Dict['brain']['rheMac2']['ac'] = ['/home/ak2267/project/EnhancerClasses/rheMac2/ac/merge_','_overlap_rh_named.bed']
PeakPath_Dict['brain']['rheMac2']['me2'] = ['/home/ak2267/project/EnhancerClasses/rheMac2/me2/merge_','_overlap_me2_rh_named.bed']
PeakPath_Dict['brain']['mm9']['ac'] = ['/home/ak2267/project/EnhancerClasses/mm9/ac/merge_','_overlap_mm_named.bed']
PeakPath_Dict['brain']['mm9']['me2'] = ['/home/ak2267/project/EnhancerClasses/mm9/me2/merge_','_overlap_me2_mm_named.bed']
PeakPath_Dict['limb']['hg19']['ac'] = ['/home/ak2267/project/EnhancerClasses/Limb/hg19/merge_','_overlap_named.bed']
PeakPath_Dict['limb']['rheMac2']['ac'] = ['/home/ak2267/project/EnhancerClasses/Limb/rheMac2/merge_','_overlap_named.bed']
PeakPath_Dict['limb']['mm9']['ac'] = ['/home/ak2267/project/EnhancerClasses/Limb/mm9/merge_','_overlap_named.bed']

# store paths to bigWig files in dictionary - all require timePoint to be inserted between 0 and 1 in array, and rep between 1 and 2
BigWigPath_Dict = {}
for tissue in ['brain','limb']:
    BigWigPath_Dict[tissue] = {}
    for species in ['hg19','rheMac2','mm9']:
        BigWigPath_Dict[tissue][species] = {}
BigWigPath_Dict['brain']['hg19']['ac'] = ['/gpfs/gibbs/pi/noonan/ak2267/NoonanBrain/hg19/ac/','_ac_rep','.bw']
BigWigPath_Dict['brain']['hg19']['me2'] = ['/gpfs/gibbs/pi/noonan/ak2267/NoonanBrain/hg19/me2/','_me2_rep','.bw']
BigWigPath_Dict['brain']['rheMac2']['ac'] = ['/gpfs/gibbs/pi/noonan/ak2267/NoonanBrain/rheMac2/ac/','_ac_rep','.bw']
BigWigPath_Dict['brain']['rheMac2']['me2'] = ['/gpfs/gibbs/pi/noonan/ak2267/NoonanBrain/rheMac2/me2/','_me2_rep','.bw']
BigWigPath_Dict['brain']['mm9']['ac'] = ['/gpfs/gibbs/pi/noonan/ak2267/NoonanBrain/mm9/ac/','_ac_rep','.bw']
BigWigPath_Dict['brain']['mm9']['me2'] = ['/gpfs/gibbs/pi/noonan/ak2267/NoonanBrain/mm9/me2/','_me2_rep','.bw']
BigWigPath_Dict['limb']['hg19']['ac'] = ['/gpfs/gibbs/pi/noonan/ak2267/NoonanLimb/hg19/','_ac_rep','.bw']
BigWigPath_Dict['limb']['rheMac2']['ac'] = ['/gpfs/gibbs/pi/noonan/ak2267/NoonanLimb/rheMac2/','_ac_rep','.bw'] # want to average e31, e32, e33 rep 1
BigWigPath_Dict['limb']['mm9']['ac'] = ['/gpfs/gibbs/pi/noonan/ak2267/NoonanLimb/mm9/','_ac_rep','.bw']

# store info on rep number for each combo
RepNumber_Dict = {}
for tissue in ['brain','limb']:
    RepNumber_Dict[tissue] = {}
RepNumber_Dict['brain']['hg19'] = [2,2,2,2]
RepNumber_Dict['brain']['rheMac2'] = [0,1,2,2]
RepNumber_Dict['brain']['mm9'] =[2,2,2,2]
RepNumber_Dict['limb']['hg19'] = [2,2,2,2]
RepNumber_Dict['limb']['rheMac2'] = [0,3,1,0] # 3 replicates are all named rep1, but for timePoints e31, e32, e33 - DOUBLE CHECK 1 REP for TP 2
RepNumber_Dict['limb']['mm9'] =[3,2,2,2]

## START PIPELINE

# store filePaths and fileEnds
workingDir = '/gpfs/gibbs/pi/noonan/ak2267/speciesPairs/peakCentric'
pathToConsensusCGIs = '/gpfs/gibbs/pi/noonan/ak2267/speciesPairs/consensusCGIs/'
# then: {speciesA}_{speciesB}/{speciesA}_{speciesB}_speciesA_reconciledCGI.bed
#    or {speciesA}_{speciesB}/{speciesA}_{speciesB}_speciesB_reconciledCGI.bed
pathToAnnotations = '/home/ak2267/genomes/RefSeq/featureAnnotations/'
# then: ${species}_allFeatures.bed
pathToFANTOM = '/home/ak2267/genomes/FANTOM/'
# then: FANTOM_TSS_${species}.bed
pathToBlacklist = '/home/ak2267/genomes/blacklist/'
# then: blacklist_${species}.bed

FilesToExclude_Dict = {}
FilesToExclude_Dict['filePath'] = {'RefSeq':'/home/ak2267/genomes/RefSeq/featureAnnotations/','FANTOM':'/home/ak2267/genomes/FANTOM/FANTOM_TSS_','blacklist':'/home/ak2267/genomes/blacklist/blacklist_'}
FilesToExclude_Dict['fileEnd'] = {'RefSeq':'_allFeatures.bed','FANTOM':'.bed','blacklist':'.bed'}

SpeciesCombos = ['hg19_rheMac2','hg19_mm9','rheMac2_mm9']

for speciesCombo in SpeciesCombos:
    for tissue in Tissues:
        for mark in Marks_Dict[tissue]:
            for timePointIndex in [0,1,2,3]:
                speciesA = speciesCombo.split('_')[0]
                speciesB = speciesCombo.split('_')[1]

                timePointA = TimePoints_Dict[tissue][speciesA][timePointIndex]
                timePointB = TimePoints_Dict[tissue][speciesB][timePointIndex]

                if timePointA != '.' and timePointB != '.':

                    outputString = []

                    #########
                    # STEP 0: append info on where to work and source files
                    label = speciesCombo +'_'+ tissue +'_'+ mark + '_' + str(timePointIndex)
                    outputString.append('cd '+workingDir+' ; source /home/ak2267/.bashrc ; source /home/ak2267/.bash_profile')
                    outputString.append('mkdir Noonan_'+label+' ; cd Noonan_'+label)

                    #print(speciesCombo+'\t'+tissue+'\t'+mark)
                    speciesA = speciesCombo.split('_')[0]
                    speciesB = speciesCombo.split('_')[1]
                    reconciledCGI_A = pathToConsensusCGIs+speciesA+'_'+speciesB+'/'+speciesA+'_'+speciesB+'_speciesA_reconciledCGI.bed'
                    reconciledCGI_B = pathToConsensusCGIs+speciesA+'_'+speciesB+'/'+speciesA+'_'+speciesB+'_speciesB_reconciledCGI.bed'
                    reconciledHuman = pathToConsensusCGIs+speciesA+'_'+speciesB+'/'+speciesA+'_'+speciesB+'_hg19_reconciledCGI.bed'
                    phastConsHuman = '/home/ak2267/genomes/phastCons/phastConsElements100way_hg19.bed'
                    ageSegmentationHuman = '/home/ak2267/project/EnhancerAge/clade_specific_elements_hg19.bed'

                    pathToPeaksA = timePointA.join(PeakPath_Dict[tissue][speciesA][mark])
                    pathToPeaksB = timePointB.join(PeakPath_Dict[tissue][speciesB][mark])
                    # these are still arrays because they need both timePoint and rep filled in
                    pathToBigWigA = BigWigPath_Dict[tissue][speciesA][mark]
                    pathToBigWigB = BigWigPath_Dict[tissue][speciesB][mark]

                    ######### (analogous to steps 1-3 in consensusCGI files)
                    # STEP 1: make files with 4 columns and CGI names (as "peak" names), with 'a' or 'b' preceding in each species
                    outputString.append('cat '+pathToPeaksA+' | awk \'{ print $1\"\\t\"$2\"\\t\"$3\"\\ta_\"NR }\' > speciesA_peaks_4col.bed')
                    outputString.append('cat '+pathToPeaksB+' | awk \'{ print $1\"\\t\"$2\"\\t\"$3\"\\tb_\"NR }\' > speciesB_peaks_4col.bed')

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

					#hg38FilesToExclude = [FilesToExclude_Dict['filePath']['RefSeq']+'hg38'+FilesToExclude_Dict['fileEnd']['RefSeq'], FilesToExclude_Dict['filePath']['FANTOM']+'hg38'+FilesToExclude_Dict['fileEnd']['FANTOM'], FilesToExclude_Dict['filePath']['blacklist']+'hg38'+FilesToExclude_Dict['fileEnd']['blacklist']]
                    hg19FilesToExclude = [FilesToExclude_Dict['filePath']['RefSeq']+'hg19'+FilesToExclude_Dict['fileEnd']['RefSeq'], FilesToExclude_Dict['filePath']['FANTOM']+'hg19'+FilesToExclude_Dict['fileEnd']['FANTOM'], FilesToExclude_Dict['filePath']['blacklist']+'hg19'+FilesToExclude_Dict['fileEnd']['blacklist']]

                    # lift both species to human

                    # remove features in both species
                    outputString.append('bedtools intersect -v -a speciesA_peaks_4col.bed -b '+' -b '.join(filesToExcludeA)+' > speciesA_peaks_4col_noFeature.bed')
                    outputString.append('bedtools intersect -v -a speciesB_peaks_4col.bed -b '+' -b '.join(filesToExcludeB)+' > speciesB_peaks_4col_noFeature.bed')

                    # store chains
                    chainAtoHg19 = '/home/ak2267/genomes/chain/'+speciesA+'ToHg19.over.chain.gz'
                    chainBtoHg19 = '/home/ak2267/genomes/chain/'+speciesB+'ToHg19.over.chain.gz'
                    chainHg19toA = '/home/ak2267/genomes/chain/hg19To'+speciesA[0].upper()+speciesA[1:]+'.over.chain.gz'
                    chainHg19toB = '/home/ak2267/genomes/chain/hg19To'+speciesB[0].upper()+speciesB[1:]+'.over.chain.gz'

                    # lift to human (if not already human)
                    if speciesA == 'hg19':
                        outputString.append('cp speciesA_peaks_4col_noFeature.bed speciesA_peaks_4col_noFeature_LOhg19.bed')
                        outputString.append('liftOver -minMatch=0.3 speciesB_peaks_4col_noFeature.bed '+chainBtoHg19+' speciesB_peaks_4col_noFeature_LOhg19.bed unMapped')
                    if speciesA != 'hg19':
                        outputString.append('liftOver -minMatch=0.3 speciesA_peaks_4col_noFeature.bed '+chainAtoHg19+' speciesA_peaks_4col_noFeature_LOhg19.bed unMapped')
                        outputString.append('liftOver -minMatch=0.3 speciesB_peaks_4col_noFeature.bed '+chainBtoHg19+' speciesB_peaks_4col_noFeature_LOhg19.bed unMapped')

                    # lift back and restrict list in human to only those that lift back to the same CGI
                    if speciesA == 'hg19':
                        outputString.append('cp speciesA_peaks_4col_noFeature_LOhg19.bed speciesA_peaks_4col_noFeature_liftedBackToA.bed')
                        outputString.append('liftOver -minMatch=0.3 speciesB_peaks_4col_noFeature_LOhg19.bed '+chainHg19toB+' speciesB_peaks_4col_noFeature_liftedBackToB.bed unMapped')
                    if speciesA != 'hg19':
                        outputString.append('liftOver -minMatch=0.3 speciesA_peaks_4col_noFeature_LOhg19.bed '+chainHg19toA+' speciesA_peaks_4col_noFeature_liftedBackToA.bed unMapped')
                        outputString.append('liftOver -minMatch=0.3 speciesB_peaks_4col_noFeature_LOhg19.bed '+chainHg19toB+' speciesB_peaks_4col_noFeature_liftedBackToB.bed unMapped')
                    outputString.append('bedtools intersect -wao -a speciesA_peaks_4col_noFeature_liftedBackToA.bed -b speciesA_peaks_4col_noFeature.bed > speciesA_checkMappingBack.txt')
                    outputString.append('bedtools intersect -wao -a speciesB_peaks_4col_noFeature_liftedBackToB.bed -b speciesB_peaks_4col_noFeature.bed > speciesB_checkMappingBack.txt')
                    outputString.append('python /gpfs/gibbs/pi/noonan/ak2267/speciesPairs/consensusCGIs/Scripts/restrictToSitesThatMapBack.py speciesA_checkMappingBack.txt > speciesA_peaks_4col_noFeature_mapBack_inAcoord.bed')
                    outputString.append('python /gpfs/gibbs/pi/noonan/ak2267/speciesPairs/consensusCGIs/Scripts/restrictToSitesThatMapBack.py speciesB_checkMappingBack.txt > speciesB_peaks_4col_noFeature_mapBack_inBcoord.bed')
                    outputString.append('python /gpfs/gibbs/pi/noonan/ak2267/speciesPairs/consensusCGIs/Scripts/restrictToLO.py speciesA_peaks_4col_noFeature_LOhg19.bed speciesA_peaks_4col_noFeature_mapBack_inAcoord.bed > speciesA_peaks_4col_noFeature_LOhg19_mapBack.bed')
                    outputString.append('python /gpfs/gibbs/pi/noonan/ak2267/speciesPairs/consensusCGIs/Scripts/restrictToLO.py speciesB_peaks_4col_noFeature_LOhg19.bed speciesB_peaks_4col_noFeature_mapBack_inBcoord.bed > speciesB_peaks_4col_noFeature_LOhg19_mapBack.bed')

                    # merge, remove features in human, sort file names to be a_# then b_#
                    outputString.append('cat speciesA_peaks_4col_noFeature_LOhg19_mapBack.bed speciesB_peaks_4col_noFeature_LOhg19_mapBack.bed > cat.bed')
                    outputString.append('sort -k1,1 -k2,2n cat.bed > sort.bed')
                    outputString.append('bedtools merge -c 4 -o collapse -i sort.bed > merge_hg19.bed')
                    outputString.append('bedtools intersect -v -a merge_hg19.bed -b '+' -b '.join(hg19FilesToExclude)+' > merge_hg19_noFeature.bed')
                    outputString.append('python /gpfs/gibbs/pi/noonan/ak2267/speciesPairs/consensusCGIs/Scripts/sortRegionNames.py merge_hg19_noFeature.bed > merge_hg19_noFeature_sortedNames.bed')

                    # lift back out to speciesA and speciesB
                    if speciesA == 'hg19':
                        outputString.append('cp merge_hg19_noFeature_sortedNames.bed mergedPeaks_liftToA.bed')
                        outputString.append('liftOver -minMatch=0.3 merge_hg19_noFeature_sortedNames.bed '+chainHg19toB+' mergedPeaks_liftToB.bed unMapped')
                    if speciesA != 'hg19':
                        outputString.append('liftOver -minMatch=0.3 merge_hg19_noFeature_sortedNames.bed '+chainHg19toA+' mergedPeaks_liftToA.bed unMapped')
                        outputString.append('liftOver -minMatch=0.3 merge_hg19_noFeature_sortedNames.bed '+chainHg19toB+' mergedPeaks_liftToB.bed unMapped')

                    # lift back to human and restrict list in speciesA and speciesB to only those that lift back to the same CGI in human (analagous to above)
                    if speciesA == 'hg19':
                        outputString.append('cp mergedPeaks_liftToA.bed mergedPeaks_liftToA_backToHg19.bed')
                        outputString.append('liftOver -minMatch=0.3 mergedPeaks_liftToB.bed '+chainBtoHg19+' mergedPeaks_liftToB_backToHg19.bed unMapped')
                    if speciesA != 'hg19':
                        outputString.append('liftOver -minMatch=0.3 mergedPeaks_liftToA.bed '+chainAtoHg19+' mergedPeaks_liftToA_backToHg19.bed unMapped')
                        outputString.append('liftOver -minMatch=0.3 mergedPeaks_liftToB.bed '+chainBtoHg19+' mergedPeaks_liftToB_backToHg19.bed unMapped')

                    outputString.append('bedtools intersect -wao -a mergedPeaks_liftToA_backToHg19.bed -b merge_hg19_noFeature_sortedNames.bed > speciesA_checkMappingBackInHg19.txt')
                    outputString.append('bedtools intersect -wao -a mergedPeaks_liftToB_backToHg19.bed -b merge_hg19_noFeature_sortedNames.bed > speciesB_checkMappingBackInHg19.txt')
                    outputString.append('python /gpfs/gibbs/pi/noonan/ak2267/speciesPairs/consensusCGIs/Scripts/restrictToSitesThatMapBack.py speciesA_checkMappingBackInHg19.txt > speciesA_merged_thatMapBackToHg19_inHg19.bed')
                    outputString.append('python /gpfs/gibbs/pi/noonan/ak2267/speciesPairs/consensusCGIs/Scripts/restrictToSitesThatMapBack.py speciesB_checkMappingBackInHg19.txt > speciesB_merged_thatMapBackToHg19_inHg19.bed')
                    outputString.append('python /gpfs/gibbs/pi/noonan/ak2267/speciesPairs/consensusCGIs/Scripts/restrictToLO.py mergedPeaks_liftToA.bed speciesA_merged_thatMapBackToHg19_inHg19.bed > merged_speciesA_mapBack.bed')
                    outputString.append('python /gpfs/gibbs/pi/noonan/ak2267/speciesPairs/consensusCGIs/Scripts/restrictToLO.py mergedPeaks_liftToB.bed speciesB_merged_thatMapBackToHg19_inHg19.bed > merged_speciesB_mapBack.bed')

                    # remove features
                    outputString.append('bedtools intersect -v -a merged_speciesA_mapBack.bed -b '+' -b '.join(filesToExcludeA)+' > mergedPeaks_liftToA_noFeature.bed')
                    outputString.append('bedtools intersect -v -a merged_speciesB_mapBack.bed -b '+' -b '.join(filesToExcludeB)+' > mergedPeaks_liftToB_noFeature.bed')

                    #########
                    # STEP 3: restrict all three lists (in A, in B, and in hg19) to those present in both A and B
                    outputString.append('python /gpfs/gibbs/pi/noonan/ak2267/speciesPairs/consensusCGIs/Scripts/restrictToLO.py mergedPeaks_liftToA_noFeature.bed mergedPeaks_liftToB_noFeature.bed > speciesA_reconciledPeaks.bed')
                    outputString.append('python /gpfs/gibbs/pi/noonan/ak2267/speciesPairs/consensusCGIs/Scripts/restrictToLO.py mergedPeaks_liftToB_noFeature.bed mergedPeaks_liftToA_noFeature.bed > speciesB_reconciledPeaks.bed')
                    outputString.append('python /gpfs/gibbs/pi/noonan/ak2267/speciesPairs/consensusCGIs/Scripts/restrictToLO.py merge_hg19_noFeature_sortedNames.bed mergedPeaks_liftToA_noFeature.bed > intermediateHg19_reconciledFile.bed')
                    outputString.append('python /gpfs/gibbs/pi/noonan/ak2267/speciesPairs/consensusCGIs/Scripts/restrictToLO.py intermediateHg19_reconciledFile.bed mergedPeaks_liftToB_noFeature.bed > human_reconciledPeaks.bed')

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


					##########
					# STEP 6:
					# run bigWigAverageOverBed to get average signal in peak intervals - THIS DIFFERS FROM ROLLER PIPELINE
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

					##########
					# STEP 7: run bedtools getfasta and faCount to get CpG numbers in reconciled CGIs, and in reconciled Peaks
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
					#	outputString.append('bedtools intersect -wao -a ' + reconciledHuman + ' -b '+pathToPhastBias+phastBias_SpeciesTranslator[speciesA]+'_allChr_phastBias_merged.bed > humanCGI_intersectPhastBias_speciesA.txt')
					#	outputString.append('bedtools intersect -wao -a ' + reconciledHuman + ' -b '+pathToPhastBias+phastBias_SpeciesTranslator[speciesB]+'_allChr_phastBias_merged.bed > humanCGI_intersectPhastBias_speciesB.txt')
					#	outputString.append('bedtools intersect -wao -a human_reconciledPeaks.bed -b '+pathToPhastBias+phastBias_SpeciesTranslator[speciesA]+'_allChr_phastBias_merged.bed > humanPeaks_intersectPhastBias_speciesA.txt')
					#	outputString.append('bedtools intersect -wao -a human_reconciledPeaks.bed -b '+pathToPhastBias+phastBias_SpeciesTranslator[speciesB]+'_allChr_phastBias_merged.bed > humanPeaks_intersectPhastBias_speciesB.txt')

					#	pythonInputFiles.append('humanCGI_intersectPhastBias_speciesA.txt')
					#	pythonInputFiles.append('humanCGI_intersectPhastBias_speciesB.txt')
					#	pythonInputFiles.append('humanPeaks_intersectPhastBias_speciesA.txt')
					#	pythonInputFiles.append('humanPeaks_intersectPhastBias_speciesB.txt')
				
					# STEP 9: python script to analyze the above intersection files & output final counts
					#echo $species $tissue $mark
                    outputString.append('python /gpfs/gibbs/pi/noonan/ak2267/speciesPairs/Scripts/countFinalOverlaps_speciesPairs_peakCentric.py '+' '.join(pythonInputFiles)+' >> '+workingDir+'/Noonan_summaryFiles/'+label+'.txt')

                    print('; '.join(outputString))
            
