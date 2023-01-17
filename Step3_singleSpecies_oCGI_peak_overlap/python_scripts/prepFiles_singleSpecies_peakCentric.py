# 10/3/22
# Purpose: make job files for running tasks from peakCentric version of singleSpecies analysis
# see bottom half of 220912_CGI_activityHeatmap.sh

Roller_speciesArray = ['rheMac10','calJac4','mm39','rn7','susScr11','canFam6','felCat9','equCab3']
Roller_tissueArray = ['brain','liver','muscle','testis']
Roller_markArray = ['H3K27ac','H3K4me3','H3K4me1']

basePath = '/gpfs/gibbs/pi/noonan/ak2267/singleSpecies/220929_activitySummary_peakCentric/'

pathToAnnotations = '/home/ak2267/genomes/RefSeq/featureAnnotations/'
# then: ${species}_allFeatures.bed
pathToPhastBias = '/gpfs/gibbs/pi/noonan/ak2267/speciesPairs/SequenceAnalysis/gBGC_analysis/phastBias_bed/'
# then: ${species}_allChr_phastBias_merged.bed - note this is in hg38 coordinates
pathToCGIs = '/home/ak2267/genomes/CGI/UCSC_AL/'
# then: ${species}_CGIsAL.bed
phastConsHuman = '/home/ak2267/genomes/phastCons/phastConsElements100way_hg38.bed'
ageSegmentationHuman = '/home/ak2267/project/EnhancerAge/clade_specific_elements_hg38.bed'
pathToRepeats = '/home/ak2267/genomes/rmsk/rmsk_'
# then: ${species}_merged.bed

hg38FeaturesArray = []
hg38FeaturesArray.append(pathToAnnotations+'hg38_allFeatures.bed')
hg38FeaturesArray.append('/home/ak2267/genomes/FANTOM/FANTOM_TSS_hg38.bed')
hg38FeaturesArray.append('/home/ak2267/genomes/blacklist/blacklist_hg38.bed')

# store usable reps for each species (i.e. excluding those with low peak numbers as previously determined)
inReps = open('/gpfs/gibbs/pi/noonan/ak2267/Roller/peaks/repsToUseRoller.txt','rt')
Rep_Dict = {}
# Rep_Dict[species_tissue_mark] = [1,2,3] # or subset if a rep is unused
for line in inReps:
    splitLine = line.strip().split()
    Rep_Dict[splitLine[0]+'_'+splitLine[1]+'_'+splitLine[2]] = splitLine[3].split(',')


for species in Roller_speciesArray:
    
    featuresArray = []
    featuresArray.append(pathToAnnotations+species+'_allFeatures.bed')
    
    if species == 'mm39':
        featuresArray.append('/home/ak2267/genomes/FANTOM/FANTOM_TSS_'+species+'.bed')
        featuresArray.append('/home/ak2267/genomes/blacklist/blacklist_'+species+'.bed')

    for tissue in Roller_tissueArray:
        for mark in Roller_markArray:
            fileStem = species + '_' + tissue + '_' + mark       
            
            chainRollerToHg38 = '/home/ak2267/genomes/chain/'+species+'ToHg38.over.chain.gz'
            chainHg38toRoller = '/home/ak2267/genomes/chain/hg38To'+species[0].upper()+species[1:]+'.over.chain.gz'
 
            # initiate arrays to collect commands and input files for python script in final step
            outputArray = []
            pythonInputFiles = []
            
            # cd, source
            outputArray.append('cd '+basePath + ' ; source /home/ak2267/.bashrc ; source /home/ak2267/.bash_profile')
    
            # STEP 1: restrict to intronic/intergenic and intersect with CGIs
    
            # restrict to intronic/intergenic in Roller species
            outputArray.append('bedtools intersect -v -a /gpfs/gibbs/pi/noonan/ak2267/Roller/peaks/intersection/'+fileStem+'_intersection.bed -b '+' -b '.join(featuresArray)+' > '+fileStem+'_noFeatures.bed')
            
            pythonInputFiles.append(fileStem+'_noFeatures.bed')
    
            # lift to human
            outputArray.append('liftOver -minMatch=0.3 '+fileStem+'_noFeatures.bed '+chainRollerToHg38+' liftToHuman/'+fileStem+'_noFeatures_LOhuman.bed unMapped')
    
            # lift back to Roller species and restrict human bed file based on overlapping original site
            outputArray.append('liftOver -minMatch=0.3 liftToHuman/'+fileStem+'_noFeatures_LOhuman.bed '+chainHg38toRoller+' liftToHuman/'+fileStem+'_noFeatures_LOhuman_mapBack.bed unMapped')
            outputArray.append('bedtools intersect -wao -a liftToHuman/'+fileStem+'_noFeatures_LOhuman_mapBack.bed -b '+fileStem+'_noFeatures.bed > liftToHuman/'+fileStem+'_checkMappingBack.txt')
            outputArray.append('python /gpfs/gibbs/pi/noonan/ak2267/speciesPairs/consensusCGIs/Scripts/restrictToSitesThatMapBack.py liftToHuman/'+fileStem+'_checkMappingBack.txt > liftToHuman/'+fileStem+'_peaksThatMapBack_inRollerCoord.bed')
            outputArray.append('python /gpfs/gibbs/pi/noonan/ak2267/speciesPairs/consensusCGIs/Scripts/restrictToLO.py liftToHuman/'+fileStem+'_noFeatures_LOhuman.bed liftToHuman/'+fileStem+'_peaksThatMapBack_inRollerCoord.bed > liftToHuman/'+fileStem+'_peaksThatMapBack_inHumanCoord.bed')
            
            pythonInputFiles.append('liftToHuman/'+fileStem+'_peaksThatMapBack_inHumanCoord.bed')
    
            # restrict to intronic/intergenic in human - this file will be used to set 1/0 column in final table
            outputArray.append('bedtools intersect -v -a liftToHuman/'+fileStem+'_peaksThatMapBack_inHumanCoord.bed -b '+' -b '.join(hg38FeaturesArray)+' > liftToHuman/'+fileStem+'_noFeaturesInHuman_inHumanCoord.bed')
            
            pythonInputFiles.append('liftToHuman/'+fileStem+'_noFeaturesInHuman_inHumanCoord.bed')
    
            # intersect with CGIs in Roller species
            # first make CGI files with name in col 4
            outputArray.append('awk \'{ print $1\"\\t\"$2\"\\t\"$3\"\\tCGI_\"NR } \' '+pathToCGIs+species+'_CGIsAL.bed > '+fileStem+'_namedCGIs.bed')
            outputArray.append('bedtools intersect -wao -a '+fileStem+'_noFeatures.bed -b '+fileStem+'_namedCGIs.bed > intersectCGIs/'+fileStem+'_intersectCGIs_RollerCoord.txt')
            
            pythonInputFiles.append('intersectCGIs/'+fileStem+'_intersectCGIs_RollerCoord.txt')
    
            # STEP 2: run featureCounts to count reads in all peaks (not just intronic/intergenic ones)
            
            #turn bed file into GTF and load Subread
            outputArray.append('python /home/ak2267/Scripts/makeGTF_nameFromCol4.py /gpfs/gibbs/pi/noonan/ak2267/Roller/peaks/intersection/'+fileStem+'_intersection.bed featureCounts/'+fileStem+'.gtf')
            outputArray.append('module load Subread/2.0.0-GCC-7.3.0-2.30')
            
            #run featureCounts and append to 'file list' for use with countRPM.py in the next step
            fileList = []
            for replicate in Rep_Dict[fileStem]:
               outputArray.append('featureCounts -t exon -g gene_id -a featureCounts/'+fileStem+'.gtf -o featureCounts/'+fileStem+'_'+str(replicate)+'.counts /gpfs/gibbs/pi/noonan/ak2267/Roller/bam/'+fileStem+'_'+str(replicate)+'_bowtie2_filtered.bam')
               fileList.append('featureCounts/'+fileStem+'_'+str(replicate)+'.counts')
               
            # generate file with RPM for each reconciledPeak (simpler than doing this in the final python script)
            outputArray.append('python /gpfs/gibbs/pi/noonan/ak2267/speciesPairs/Scripts/countRPM_general.py '+str(len(fileList))+' '+' '.join(fileList)+' > featureCounts/'+fileStem+'.signal')

            pythonInputFiles.append('featureCounts/'+fileStem+'.signal')

            # STEP 3: count CpGs in peaks AND CGIs using bedtools getFasta and faCount
            outputArray.append('bedtools getfasta -name -fi /gpfs/gibbs/pi/noonan/ak2267/genomes/'+species+'.fa -bed '+fileStem+'_noFeatures.bed > faCount/'+fileStem+'_Peaks.fa ; faCount faCount/'+fileStem+'_Peaks.fa > faCount/'+fileStem+'_Peaks.faCount')
            outputArray.append('bedtools getfasta -name -fi /gpfs/gibbs/pi/noonan/ak2267/genomes/'+species+'.fa -bed '+fileStem+'_namedCGIs.bed > faCount/'+fileStem+'_CGIs.fa ; faCount faCount/'+fileStem+'_CGIs.fa > faCount/'+fileStem+'_CGIs.faCount')
            pythonInputFiles.append('faCount/'+fileStem+'_Peaks.faCount')
            pythonInputFiles.append('faCount/'+fileStem+'_CGIs.faCount')

            # STEP 4: intersect human lifted peaks AND CGIs with age segmentation
            outputArray.append('bedtools intersect -wao -a liftToHuman/'+fileStem+'_noFeaturesInHuman_inHumanCoord.bed -b '+ageSegmentationHuman+' > intersectAge/'+fileStem+'_Peaks_intersectAgeSegmentation.txt')
            # for CGIs, requires lifting CGIs to human first (and going through process of making sure they map back)
            outputArray.append('liftOver -minMatch=0.3 '+fileStem+'_namedCGIs.bed '+chainRollerToHg38+' intersectAge/'+fileStem+'_CGIsAL_liftHuman.bed unMapped')
            outputArray.append('liftOver -minMatch=0.3 intersectAge/'+fileStem+'_CGIsAL_liftHuman.bed '+chainHg38toRoller+' intersectAge/'+fileStem+'_CGIsAL_liftHuman_mapBack.bed unMapped')
            outputArray.append('bedtools intersect -wao -a intersectAge/'+fileStem+'_CGIsAL_liftHuman_mapBack.bed -b '+fileStem+'_namedCGIs.bed > intersectAge/'+fileStem+'_checkMappingBack.txt')
            outputArray.append('python /gpfs/gibbs/pi/noonan/ak2267/speciesPairs/consensusCGIs/Scripts/restrictToSitesThatMapBack.py intersectAge/'+fileStem+'_checkMappingBack.txt > intersectAge/'+fileStem+'_CGIsThatMapBack_inRollerCoord.bed')
            outputArray.append('python /gpfs/gibbs/pi/noonan/ak2267/speciesPairs/consensusCGIs/Scripts/restrictToLO.py intersectAge/'+fileStem+'_CGIsAL_liftHuman.bed intersectAge/'+fileStem+'_CGIsThatMapBack_inRollerCoord.bed > intersectAge/'+fileStem+'_CGIsThatMapBack_inHumanCoord.bed')
            outputArray.append('bedtools intersect -wao -a intersectAge/'+fileStem+'_CGIsThatMapBack_inHumanCoord.bed -b '+ageSegmentationHuman+' > intersectAge/'+fileStem+'_CGIs_intersectAgeSegmentation.txt')
            
            pythonInputFiles.append('intersectAge/'+fileStem+'_Peaks_intersectAgeSegmentation.txt')
            pythonInputFiles.append('intersectAge/'+fileStem+'_CGIs_intersectAgeSegmentation.txt')
            
            # STEP 5: intersect human lifted peaks AND CGIs with phastCons
            outputArray.append('bedtools intersect -wao -a liftToHuman/'+fileStem+'_noFeaturesInHuman_inHumanCoord.bed -b '+phastConsHuman+' > intersectPhastCons/'+fileStem+'_Peaks_intersectPhastCons.txt')
            pythonInputFiles.append('intersectPhastCons/'+fileStem+'_Peaks_intersectPhastCons.txt')
            outputArray.append('bedtools intersect -wao -a intersectAge/'+fileStem+'_CGIsThatMapBack_inHumanCoord.bed -b '+phastConsHuman+' > intersectPhastCons/'+fileStem+'_CGIs_intersectPhastCons.txt')
            pythonInputFiles.append('intersectPhastCons/'+fileStem+'_CGIs_intersectPhastCons.txt')
                        
            # STEP 6: intersect with repeats
            outputArray.append('bedtools intersect -wao -a '+fileStem+'_noFeatures.bed -b '+pathToRepeats+species+'_merged.bed > intersectRepeats/'+fileStem+'_intersectRepeats.txt')
            pythonInputFiles.append('intersectRepeats/'+fileStem+'_intersectRepeats.txt')
            
            # STEP 7: integrate all info with python script
            outputArray.append('python /gpfs/gibbs/pi/noonan/ak2267/singleSpecies/Scripts/makeSummary_singleSpecies_peakCentric.py '+' '.join(pythonInputFiles)+' > summaryFiles_singleSpeciesPeakCentric/'+fileStem+'.txt')

            # WRITE TO OUTPUT
            print(' ; '.join(outputArray))
            
# NOONAN   
            
Noonan_speciesArray = ['hg19','rheMac2','mm9']
Noonan_tissueDict = {'brain':['ac','me2'], 'limb':['ac']}

RepNumber_Dict = {}
for tissue in ['brain','limb']:
    RepNumber_Dict[tissue] = {}
RepNumber_Dict['brain']['hg19'] = [2,2,2,2]
RepNumber_Dict['brain']['rheMac2'] = [0,1,2,2]
RepNumber_Dict['brain']['mm9'] =[2,2,2,2]
RepNumber_Dict['limb']['hg19'] = [2,2,2,2]
RepNumber_Dict['limb']['rheMac2'] = [0,3,1,0] # 3 replicates are all named rep1, but for timePoints e31, e32, e33 - DOUBLE CHECK 1 REP for TP 2
RepNumber_Dict['limb']['mm9'] =[3,2,2,2]

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


phastConsHuman = '/home/ak2267/genomes/phastCons/phastConsElements100way_hg19.bed'
ageSegmentationHuman = '/home/ak2267/project/EnhancerAge/clade_specific_elements_hg19.bed'

hg19FeaturesArray = []
hg19FeaturesArray.append(pathToAnnotations+'hg19_allFeatures.bed')
hg19FeaturesArray.append('/home/ak2267/genomes/FANTOM/FANTOM_TSS_hg19.bed')
hg19FeaturesArray.append('/home/ak2267/genomes/blacklist/blacklist_hg19.bed')
            
for species in Noonan_speciesArray:

    featuresArray = []
    featuresArray.append(pathToAnnotations+species+'_allFeatures.bed')
    
    if species == 'mm9' or species == 'hg19':
        featuresArray.append('/home/ak2267/genomes/FANTOM/FANTOM_TSS_'+species+'.bed')
        featuresArray.append('/home/ak2267/genomes/blacklist/blacklist_'+species+'.bed')

    for tissue in Noonan_tissueDict:
        for mark in Noonan_tissueDict[tissue]:
            for timePointIndex in [0, 1, 2, 3]:
                if RepNumber_Dict[tissue][species][timePointIndex] > 0:
                
                    timePoint = TimePoints_Dict[tissue][species][timePointIndex]
                    
                    fileStem = species + '_' + tissue + '_' + mark + '_' + str(timePointIndex)
        
                    chainNoonanToHg19 = '/home/ak2267/genomes/chain/'+species+'ToHg19.over.chain.gz'
                    chainHg19toNoonan = '/home/ak2267/genomes/chain/hg19To'+species[0].upper()+species[1:]+'.over.chain.gz'

                    # initiate arrays to collect commands and input files for python script in final step
                    outputArray = []
                    pythonInputFiles = []
        
                    # cd, source
                    outputArray.append('cd '+basePath + ' ; source /home/ak2267/.bashrc ; source /home/ak2267/.bash_profile')

                    # STEP 1: restrict to intronic/intergenic and intersect with CGIs

                    # restrict to intronic/intergenic in Noonan species
                    peakBedFile = timePoint.join(PeakPath_Dict[tissue][species][mark])
                    outputArray.append('bedtools intersect -v -a '+ peakBedFile +' -b '+' -b '.join(featuresArray)+' > '+fileStem+'_noFeatures.bed')
        
                    pythonInputFiles.append(fileStem+'_noFeatures.bed')

                    # lift to human
                    if species != 'hg19':
                        outputArray.append('liftOver -minMatch=0.3 '+fileStem+'_noFeatures.bed '+chainNoonanToHg19+' liftToHuman/'+fileStem+'_noFeatures_LOhuman.bed unMapped')
                    elif species == 'hg19':
                        outputArray.append('cp '+fileStem+'_noFeatures.bed liftToHuman/'+fileStem+'_noFeatures_LOhuman.bed')

                    # lift back to Noonan species and restrict human bed file based on overlapping original site
                    if species != 'hg19':
                        outputArray.append('liftOver -minMatch=0.3 liftToHuman/'+fileStem+'_noFeatures_LOhuman.bed '+chainHg19toNoonan+' liftToHuman/'+fileStem+'_noFeatures_LOhuman_mapBack.bed unMapped')
                        outputArray.append('bedtools intersect -wao -a liftToHuman/'+fileStem+'_noFeatures_LOhuman_mapBack.bed -b '+fileStem+'_noFeatures.bed > liftToHuman/'+fileStem+'_checkMappingBack.txt')
                        outputArray.append('python /gpfs/gibbs/pi/noonan/ak2267/speciesPairs/consensusCGIs/Scripts/restrictToSitesThatMapBack.py liftToHuman/'+fileStem+'_checkMappingBack.txt > liftToHuman/'+fileStem+'_peaksThatMapBack_inNoonanCoord.bed')
                        outputArray.append('python /gpfs/gibbs/pi/noonan/ak2267/speciesPairs/consensusCGIs/Scripts/restrictToLO.py liftToHuman/'+fileStem+'_noFeatures_LOhuman.bed liftToHuman/'+fileStem+'_peaksThatMapBack_inNoonanCoord.bed > liftToHuman/'+fileStem+'_peaksThatMapBack_inHumanCoord.bed')
                    elif species == 'hg19':
                        outputArray.append('cp '+fileStem+'_noFeatures.bed liftToHuman/'+fileStem+'_peaksThatMapBack_inHumanCoord.bed')
        
                    pythonInputFiles.append('liftToHuman/'+fileStem+'_peaksThatMapBack_inHumanCoord.bed')

                    # restrict to intronic/intergenic in human - this file will be used to set 1/0 column in final table
                    if species != 'hg19':
                        outputArray.append('bedtools intersect -v -a liftToHuman/'+fileStem+'_peaksThatMapBack_inHumanCoord.bed -b '+' -b '.join(hg19FeaturesArray)+' > liftToHuman/'+fileStem+'_noFeaturesInHuman_inHumanCoord.bed')
                    elif species == 'hg19':
                        outputArray.append('cp liftToHuman/'+fileStem+'_peaksThatMapBack_inHumanCoord.bed liftToHuman/'+fileStem+'_noFeaturesInHuman_inHumanCoord.bed')
        
                    pythonInputFiles.append('liftToHuman/'+fileStem+'_noFeaturesInHuman_inHumanCoord.bed')

                    # intersect with CGIs in Noonan species
                    # first make CGI files with name in col 4
                    outputArray.append('awk \'{ print $1\"\\t\"$2\"\\t\"$3\"\\tCGI_\"NR } \' '+pathToCGIs+species+'_CGIsAL.bed > '+fileStem+'_namedCGIs.bed')
                    outputArray.append('bedtools intersect -wao -a '+fileStem+'_noFeatures.bed -b '+fileStem+'_namedCGIs.bed > intersectCGIs/'+fileStem+'_intersectCGIs_NoonanCoord.txt')
        
                    pythonInputFiles.append('intersectCGIs/'+fileStem+'_intersectCGIs_NoonanCoord.txt')

                    # STEP 2: run bigWigAverageOverBed to count SIGNAL in all peaks (not just intronic/intergenic ones)
        
                    # run bigWigAverageOverBed and append to 'file list' for use with countBigWigSignal.py in the next step
                    
                    # change mouse timePoint name (bigWigs are named differently than peak files)
                    timePointBigWig = timePoint
                    fixMouseTimePoints = {'17F':'e17F','17O':'e17O'}
                    if timePoint in fixMouseTimePoints:
                    	timePointBigWig = fixMouseTimePoints[timePoint]
                    pathToBigWig = BigWigPath_Dict[tissue][species][mark]
                    
                    # command
                    fileList = []
                    if timePointIndex == 1 and tissue == 'limb' and species == 'rheMac2':
                        outputArray.append('bigWigAverageOverBed '+pathToBigWig[0]+'e31'+pathToBigWig[1]+'1'+pathToBigWig[2]+' '+peakBedFile+' '+fileStem+'_1.counts')
                        outputArray.append('bigWigAverageOverBed '+pathToBigWig[0]+'e32'+pathToBigWig[1]+'1'+pathToBigWig[2]+' '+peakBedFile+' '+fileStem+'_2.counts')
                        outputArray.append('bigWigAverageOverBed '+pathToBigWig[0]+'e33'+pathToBigWig[1]+'1'+pathToBigWig[2]+' '+peakBedFile+' '+fileStem+'_3.counts')
                        fileList.append(fileStem+'_1.counts')
                        fileList.append(fileStem+'_2.counts')
                        fileList.append(fileStem+'_3.counts')

                    else:
                        for replicate in range(1,RepNumber_Dict[tissue][species][timePointIndex]+1):
                            outputArray.append('bigWigAverageOverBed '+pathToBigWig[0]+timePointBigWig+pathToBigWig[1]+str(replicate)+pathToBigWig[2]+' '+peakBedFile+' '+fileStem+'_'+str(replicate)+'.counts')
                            fileList.append(fileStem+'_'+str(replicate)+'.counts')
           
                    # generate file with RPM for each reconciledPeak (simpler than doing this in the final python script)
                    outputArray.append('python /gpfs/gibbs/pi/noonan/ak2267/speciesPairs/Scripts/countBigWigSignal_general.py '+str(len(fileList))+' '+' '.join(fileList)+' > featureCounts/'+fileStem+'.signal')

                    pythonInputFiles.append('featureCounts/'+fileStem+'.signal')

                    # STEP 3: count CpGs in peaks AND CGIs using bedtools getFasta and faCount
                    outputArray.append('bedtools getfasta -name -fi /gpfs/gibbs/pi/noonan/ak2267/genomes/'+species+'.fa -bed '+fileStem+'_noFeatures.bed > faCount/'+fileStem+'_Peaks.fa ; faCount faCount/'+fileStem+'_Peaks.fa > faCount/'+fileStem+'_Peaks.faCount')
                    outputArray.append('bedtools getfasta -name -fi /gpfs/gibbs/pi/noonan/ak2267/genomes/'+species+'.fa -bed '+fileStem+'_namedCGIs.bed > faCount/'+fileStem+'_CGIs.fa ; faCount faCount/'+fileStem+'_CGIs.fa > faCount/'+fileStem+'_CGIs.faCount')
                    pythonInputFiles.append('faCount/'+fileStem+'_Peaks.faCount')
                    pythonInputFiles.append('faCount/'+fileStem+'_CGIs.faCount')

                    # STEP 4: intersect human lifted peaks AND CGIs with age segmentation
                    if species != 'hg19':
                        # for CGIs, requires lifting CGIs to human first (and going through process of making sure they map back)
                        outputArray.append('liftOver -minMatch=0.3 '+fileStem+'_namedCGIs.bed '+chainNoonanToHg19+' intersectAge/'+fileStem+'_CGIsAL_liftHuman.bed unMapped')
                        outputArray.append('liftOver -minMatch=0.3 intersectAge/'+fileStem+'_CGIsAL_liftHuman.bed '+chainHg19toNoonan+' intersectAge/'+fileStem+'_CGIsAL_liftHuman_mapBack.bed unMapped')
                        outputArray.append('bedtools intersect -wao -a intersectAge/'+fileStem+'_CGIsAL_liftHuman_mapBack.bed -b '+fileStem+'_namedCGIs.bed > intersectAge/'+fileStem+'_checkMappingBack.txt')
                        outputArray.append('python /gpfs/gibbs/pi/noonan/ak2267/speciesPairs/consensusCGIs/Scripts/restrictToSitesThatMapBack.py intersectAge/'+fileStem+'_checkMappingBack.txt > intersectAge/'+fileStem+'_CGIsThatMapBack_inNoonanCoord.bed')
                        outputArray.append('python /gpfs/gibbs/pi/noonan/ak2267/speciesPairs/consensusCGIs/Scripts/restrictToLO.py intersectAge/'+fileStem+'_CGIsAL_liftHuman.bed intersectAge/'+fileStem+'_CGIsThatMapBack_inNoonanCoord.bed > intersectAge/'+fileStem+'_CGIsThatMapBack_inHumanCoord.bed')
                    elif species == 'hg19':
                        outputArray.append('cp '+fileStem+'_namedCGIs.bed intersectAge/'+fileStem+'_CGIsThatMapBack_inHumanCoord.bed')
                    
                    outputArray.append('bedtools intersect -wao -a liftToHuman/'+fileStem+'_noFeaturesInHuman_inHumanCoord.bed -b '+ageSegmentationHuman+' > intersectAge/'+fileStem+'_Peaks_intersectAgeSegmentation.txt')
                    outputArray.append('bedtools intersect -wao -a intersectAge/'+fileStem+'_CGIsThatMapBack_inHumanCoord.bed -b '+ageSegmentationHuman+' > intersectAge/'+fileStem+'_CGIs_intersectAgeSegmentation.txt')
        
                    pythonInputFiles.append('intersectAge/'+fileStem+'_Peaks_intersectAgeSegmentation.txt')
                    pythonInputFiles.append('intersectAge/'+fileStem+'_CGIs_intersectAgeSegmentation.txt')
        
                    # STEP 5: intersect human lifted peaks AND CGIs with phastCons
                    outputArray.append('bedtools intersect -wao -a liftToHuman/'+fileStem+'_noFeaturesInHuman_inHumanCoord.bed -b '+phastConsHuman+' > intersectPhastCons/'+fileStem+'_Peaks_intersectPhastCons.txt')
                    pythonInputFiles.append('intersectPhastCons/'+fileStem+'_Peaks_intersectPhastCons.txt')
                    outputArray.append('bedtools intersect -wao -a intersectAge/'+fileStem+'_CGIsThatMapBack_inHumanCoord.bed -b '+phastConsHuman+' > intersectPhastCons/'+fileStem+'_CGIs_intersectPhastCons.txt')
                    pythonInputFiles.append('intersectPhastCons/'+fileStem+'_CGIs_intersectPhastCons.txt')
                    
                    # STEP 6: intersect with repeats
                    outputArray.append('bedtools intersect -wao -a '+fileStem+'_noFeatures.bed -b '+pathToRepeats+species+'_merged.bed > intersectRepeats/'+fileStem+'_intersectRepeats.txt')
                    pythonInputFiles.append('intersectRepeats/'+fileStem+'_intersectRepeats.txt')
        
                    # STEP 7: integrate all info with python script
                    outputArray.append('python /gpfs/gibbs/pi/noonan/ak2267/singleSpecies/Scripts/makeSummary_singleSpecies_peakCentric.py '+' '.join(pythonInputFiles)+' > summaryFiles_singleSpeciesPeakCentric/'+fileStem+'.txt')

                    # WRITE TO OUTPUT
                    print(' ; '.join(outputArray))

