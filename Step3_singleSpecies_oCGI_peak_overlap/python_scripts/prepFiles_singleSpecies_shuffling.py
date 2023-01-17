# 1/12/23
# Purpose: make job files for reshuffling test of enrichment of peaks on oCGIs

#########################
###### CGI-CENTRIC ######
#########################

# Job file 1: generate oCGIs in each species that reciprocally lift to human
# Job file 2:
# reshuffle onto same genome background
# lift to hg38 and back, only keeping ones that lift
# count sites with peak in each tissue x mark and write to separate output file for each tissue x mark
# ALSO lift sites to rheMac2 or mm9 as necessary in order to overlap with marks in those genomes

# Then return to 220912_CGI_singleSpecies_activity.sh to summarize using summarizeReshuffling.py

RollerSpecies = ['rheMac10', 'calJac4', 'mm39', 'rn7', 'susScr11', 'canFam6', 'felCat9', 'equCab3']
NoonanSpecies = ['hg19', 'rheMac2', 'mm9']

outJobFile1 = open('230112_step1_get_oCGIs_jobFile.txt', 'wt')
outJobFile2 = open('230112_step2_do_shuffling_jobFile.txt', 'wt')

# JOB FILE 1:

# I already made bed files of oCGIs that lift reciprocally to hg38, in 220912_CGI_singleSpecies_activity.sh
# see section STEP 0 in CGI-centric pipeline and STEP 4 where I intersect with phastCons elements (which requires lifting to hg38 first)
# /gpfs/gibbs/pi/noonan/ak2267/singleSpecies/220919_activitySummary/intersectPhastCons/${species}_CGIsAL_noFeatures_LOhg38_uniqueMapping_noHg38features.bed
# but this is in hg38, so need to run restrictToLO.py to get list in each original genome, by restricting this bed file:
# /gpfs/gibbs/pi/noonan/ak2267/singleSpecies/220919_activitySummary/${species}_CGIsAL_noFeatures.bed

for species in RollerSpecies:
    
    outJobFile1.write('python /gpfs/gibbs/pi/noonan/ak2267/speciesPairs/consensusCGIs/Scripts/restrictToLO.py /gpfs/gibbs/pi/noonan/ak2267/singleSpecies/220919_activitySummary/'+species+'_CGIsAL_noFeatures.bed /gpfs/gibbs/pi/noonan/ak2267/singleSpecies/220919_activitySummary/intersectPhastCons/'+species+'_CGIsAL_noFeatures_LOhg38_uniqueMapping_noHg38features.bed > '+species+'_oCGIs_noFeatures.bed'+'\n')

# I also lifted from hg38 / rheMac10 / mm39 to hg19 / rheMac2 / mm9 - in bed file:
# /gpfs/gibbs/pi/noonan/ak2267/singleSpecies/220919_activitySummary/${species}_CGIsAL_noFeatures.bed
# but not that this has NOT been restricted based on not overlapping hg38 features - will need to do here

for species in NoonanSpecies:

    if species == 'hg19':
        outJobFile1.write('cp /gpfs/gibbs/pi/noonan/ak2267/singleSpecies/220919_activitySummary/'+species+'_CGIsAL_noFeatures.bed hg19_oCGIs_noFeatures.bed'+'\n')

    else:
        if species == 'rheMac2':
            species_otherVersion = 'rheMac10'
        
        if species == 'mm9':
            species_otherVersion = 'mm39'

        outJobFile1.write('python /gpfs/gibbs/pi/noonan/ak2267/speciesPairs/consensusCGIs/Scripts/restrictToLO.py /gpfs/gibbs/pi/noonan/ak2267/singleSpecies/220919_activitySummary/'+species+'_CGIsAL_noFeatures.bed /gpfs/gibbs/pi/noonan/ak2267/singleSpecies/220919_activitySummary/intersectPhastCons/'+species_otherVersion+'_CGIsAL_noFeatures_LOhg38_uniqueMapping_noHg38features.bed > '+species+'_oCGIs_noFeatures.bed'+'\n')

outJobFile1.close()


# JOB FILE 2: shuffling - these are the steps:
# 1) bedtools shuffle the oCGIs in each Roller species - exclude intervals in ${species}_merge.bed
# 2) lift these new intervals to hg38 and back, keeping only those that reciprocally lift AND do not have feature in hg38
# 3) for each tissue x mark, intersect this shuffle set restricted as described in step 2 with peaks
# 4) write number (numerator and denominator) to countsWithPeak/${species}_${tissue}_${mark}_countsWithPeak.txt

# special additional steps for Noonan genomes (rheMac10/rheMac2, mm39/mm9) lift the shuffle set in hg38 coordinates to hg19
# 5) lift shuffled set rheMac10 --> rheMac8 --> rheMac2, and mm39 --> mm10 --> mm9
# 6) lift back and restrict to reciprocal mapping sites
# 7) intersect with Noonan peaks for each tissue x mark
# 8) write number (numerator and denominator) to countsWithPeak/{species}_${tissue}_${mark}_countsWithPeak.txt

# for hg19: simply shuffle within hg19

# For all species, do 20000 shuffles but split them into x separate jobs writing into y separate intermediate files (x * y = 20000)

for species in RollerSpecies:

    # set number of jobs to write, and number of permutations per job
    totalPermutations = 20000
    
    if species in ['susScr11', 'felCat9']:
        numberOfJobs = 200
        
    elif species in ['canFam6', 'equCab3']:
        numberOfJobs = 100
        
    else:
        numberOfJobs = 50
        
    permutationsPerJob = totalPermutations / numberOfJobs

    # do this repeatedly for numberOfJobs times - use i to keep intermediate files distinct
    for i in range(1,numberOfJobs+1):
    
        outArray = []
    
        # cd and source
        outArray.append('cd /gpfs/gibbs/pi/noonan/ak2267/singleSpecies/shuffling ; source /home/ak2267/.bashrc ; source /home/ak2267/.bash_profile')
            
        # shuffle oCGIs
        outArray.append('for i in {1..'+str(permutationsPerJob)+'}; do bedtools shuffle -excl '+species+'_merge.bed -chrom -noOverlapping -i '+species+'_oCGIs_noFeatures.bed -g /home/ak2267/genomes/chrom.sizes/'+species+'.chrom.sizes > intermediate_'+species+'/'+species+'_shuffled_jobNum'+str(i)+'.bed')
        # note that I'm breaking here in the middle of a unix for loop
        
        # remove previously used file for collecting oCGIs with ANY peak to start fresh
        outArray.append('rm intermediate_'+species+'/'+species+'_jobNum'+str(i)+'_allIntersects.bed')
        
        # lift to hg38 and back, with step in hg38 where I restrict on hg38 features
        outArray.append('liftOver -minMatch=0.3 intermediate_'+species+'/'+species+'_shuffled_jobNum'+str(i)+'.bed /home/ak2267/genomes/chain/'+species+'ToHg38.over.chain.gz intermediate_'+species+'/'+species+'_shuffled_jobNum'+str(i)+'_LOhg38.bed unMapped')
        outArray.append('bedtools intersect -v -a intermediate_'+species+'/'+species+'_shuffled_jobNum'+str(i)+'_LOhg38.bed -b hg38_merge.bed > intermediate_'+species+'/'+species+'_shuffled_jobNum'+str(i)+'_LOhg38_noFeatures.bed')
        outArray.append('liftOver -minMatch=0.3 intermediate_'+species+'/'+species+'_shuffled_jobNum'+str(i)+'_LOhg38_noFeatures.bed /home/ak2267/genomes/chain/hg38To'+species[0].upper()+species[1:]+'.over.chain.gz intermediate_'+species+'/'+species+'_shuffled_jobNum'+str(i)+'_LOhg38_mapBack.bed unMapped')
        
        # filter list in RollerSpecies genome to those that lift reciprocally
        outArray.append('bedtools intersect -wao -a intermediate_'+species+'/'+species+'_shuffled_jobNum'+str(i)+'_LOhg38_mapBack.bed -b intermediate_'+species+'/'+species+'_shuffled_jobNum'+str(i)+'.bed > intermediate_'+species+'/'+species+'_shuffled_jobNum'+str(i)+'_checkMapBack.txt')
        outArray.append('python /gpfs/gibbs/pi/noonan/ak2267/speciesPairs/consensusCGIs/Scripts/restrictToSitesThatMapBack.py intermediate_'+species+'/'+species+'_shuffled_jobNum'+str(i)+'_checkMapBack.txt > intermediate_'+species+'/'+species+'_shuffled_jobNum'+str(i)+'_mapBothWays.bed')
        
        # get denominator of sites that lifted reciprocally
        outArray.append('denominator=$(cat intermediate_'+species+'/'+species+'_shuffled_jobNum'+str(i)+'_mapBothWays.bed | wc -l)')
        
        # for each tissue x mark, intersect with peaks, and write numbers to ${species}_${tissue}_${mark}_countsWithPeak.txt
        for tissue in ['brain', 'liver', 'muscle', 'testis']:
            for mark in ['H3K4me3', 'H3K27ac', 'H3K4me1']:
                outArray.append('withPeak=$(bedtools intersect -wa -u -a intermediate_'+species+'/'+species+'_shuffled_jobNum'+str(i)+'_mapBothWays.bed -b /gpfs/gibbs/pi/noonan/ak2267/Roller/peaks/intersection/'+species+'_'+tissue+'_'+mark+'_intersection.bed | wc -l)')
                outArray.append('echo -e $withPeak\"\\t\"$denominator >> countsWithPeak/'+species+'_'+tissue+'_'+mark+'_countsWithPeak_jobNum'+str(i)+'.txt')
                
                # write to files for ANY mark calculations
                outArray.append('bedtools intersect -wa -a intermediate_'+species+'/'+species+'_shuffled_jobNum'+str(i)+'_mapBothWays.bed -b /gpfs/gibbs/pi/noonan/ak2267/Roller/peaks/intersection/'+species+'_'+tissue+'_'+mark+'_intersection.bed >> intermediate_'+species+'/'+species+'_jobNum'+str(i)+'_allIntersects.bed')
        
        # if species also has Noonan data, do steps 5-8
        if species not in ['rheMac10', 'mm39']:
        
            # count number with ANY mark
            outArray.append('allPeaks=$(cut -f 4 intermediate_'+species+'/'+species+'_jobNum'+str(i)+'_allIntersects.bed | sort -u | wc -l)')
            outArray.append('echo -e $allPeaks\"\\t\"$denominator >> countsWithPeak/'+species+'_all_all_countsWithPeak_jobNum'+str(i)+'.txt')
        
            # close for loop
            outArray.append('done')
        
        else:
            
            # lift to rheMac8 --> rheMac2 or mm10 --> mm9
            Lift1_Dict = {'rheMac10':'rheMac8', 'mm39':'mm10'}
            Lift2_Dict = {'rheMac10':'rheMac2', 'mm39':'mm9'}
            
            outArray.append('liftOver -minMatch=0.8 intermediate_'+species+'/'+species+'_shuffled_jobNum'+str(i)+'_mapBothWays.bed /home/ak2267/genomes/chain/'+species+'To'+Lift1_Dict[species][0].upper()+Lift1_Dict[species][1:]+'.over.chain.gz intermediate_'+species+'/'+species+'_shuffled_jobNum'+str(i)+'_mapBothWays_Lift1ToNoonan.bed unMapped')
            outArray.append('liftOver -minMatch=0.8 intermediate_'+species+'/'+species+'_shuffled_jobNum'+str(i)+'_mapBothWays_Lift1ToNoonan.bed /home/ak2267/genomes/chain/'+Lift1_Dict[species]+'To'+Lift2_Dict[species][0].upper()+Lift2_Dict[species][1:]+'.over.chain.gz intermediate_'+species+'/'+species+'_shuffled_jobNum'+str(i)+'_mapBothWays_Lift2ToNoonan.bed unMapped')
            
            # lift back to rheMac10 / mm39
            outArray.append('liftOver -minMatch=0.8 intermediate_'+species+'/'+species+'_shuffled_jobNum'+str(i)+'_mapBothWays_Lift2ToNoonan.bed /home/ak2267/genomes/chain/'+Lift2_Dict[species]+'To'+Lift1_Dict[species][0].upper()+Lift1_Dict[species][1:]+'.over.chain.gz intermediate_'+species+'/'+species+'_shuffled_jobNum'+str(i)+'_mapBothWays_Lift1FromNoonan.bed unMapped')
            outArray.append('liftOver -minMatch=0.8 intermediate_'+species+'/'+species+'_shuffled_jobNum'+str(i)+'_mapBothWays_Lift1FromNoonan.bed /home/ak2267/genomes/chain/'+Lift1_Dict[species]+'To'+species[0].upper()+species[1:]+'.over.chain.gz intermediate_'+species+'/'+species+'_shuffled_jobNum'+str(i)+'_mapBothWays_Lift2FromNoonan.bed unMapped')
            
            # restrict to reciprocally mapping sites
            outArray.append('bedtools intersect -wao -a intermediate_'+species+'/'+species+'_shuffled_jobNum'+str(i)+'_mapBothWays_Lift2FromNoonan.bed -b intermediate_'+species+'/'+species+'_shuffled_jobNum'+str(i)+'_mapBothWays.bed > intermediate_'+species+'/'+species+'_shuffled_jobNum'+str(i)+'_checkMappingToNoonan.txt')
            outArray.append('python /gpfs/gibbs/pi/noonan/ak2267/speciesPairs/consensusCGIs/Scripts/restrictToSitesThatMapBack.py intermediate_'+species+'/'+species+'_shuffled_jobNum'+str(i)+'_checkMappingToNoonan.txt > intermediate_'+species+'/'+species+'_shuffled_jobNum'+str(i)+'_mapToNoonan_inRollerGenome.bed')
            outArray.append('python /gpfs/gibbs/pi/noonan/ak2267/speciesPairs/consensusCGIs/Scripts/restrictToLO.py intermediate_'+species+'/'+species+'_shuffled_jobNum'+str(i)+'_mapBothWays_Lift2ToNoonan.bed intermediate_'+species+'/'+species+'_shuffled_jobNum'+str(i)+'_mapToNoonan_inRollerGenome.bed > intermediate_'+species+'/'+species+'_shuffled_jobNum'+str(i)+'_mapToNoonan_inNoonanGenome.bed')
            
            # for each tissue x mark, intersect with peaks, and write numbers to ${species}_${tissue}_${mark}_countsWithPeak.txt
            outArray.append('withPeak=$(bedtools intersect -wa -u -a intermediate_'+species+'/'+species+'_shuffled_jobNum'+str(i)+'_mapToNoonan_inNoonanGenome.bed -b /home/ak2267/project/EnhancerClasses/'+Lift2_Dict[species]+'/ac/merge_*_overlap*_named.bed | wc -l)')
            outArray.append('echo -e $withPeak\"\\t\"$denominator >> countsWithPeak/'+species+'_devBrain_H3K27ac_countsWithPeak_jobNum'+str(i)+'.txt')
            outArray.append('withPeak=$(bedtools intersect -wa -u -a intermediate_'+species+'/'+species+'_shuffled_jobNum'+str(i)+'_mapToNoonan_inNoonanGenome.bed -b /home/ak2267/project/EnhancerClasses/'+Lift2_Dict[species]+'/me2/merge_*_overlap*_named.bed | wc -l)')
            outArray.append('echo -e $withPeak\"\\t\"$denominator >> countsWithPeak/'+species+'_devBrain_H3K4me2_countsWithPeak_jobNum'+str(i)+'.txt')
            outArray.append('withPeak=$(bedtools intersect -wa -u -a intermediate_'+species+'/'+species+'_shuffled_jobNum'+str(i)+'_mapToNoonan_inNoonanGenome.bed -b /home/ak2267/project/EnhancerClasses/Limb/'+Lift2_Dict[species]+'/merge_*_overlap_named.bed | wc -l)')
            outArray.append('echo -e $withPeak\"\\t\"$denominator >> countsWithPeak/'+species+'_devLimb_H3K27ac_countsWithPeak_jobNum'+str(i)+'.txt')
            
            # write to files for ANY mark calculations
            outArray.append('bedtools intersect -wa -a intermediate_'+species+'/'+species+'_shuffled_jobNum'+str(i)+'_mapToNoonan_inNoonanGenome.bed -b /home/ak2267/project/EnhancerClasses/'+Lift2_Dict[species]+'/ac/merge_*_overlap*_named.bed >> intermediate_'+species+'/'+species+'_jobNum'+str(i)+'_allIntersects.bed')
            outArray.append('bedtools intersect -wa -a intermediate_'+species+'/'+species+'_shuffled_jobNum'+str(i)+'_mapToNoonan_inNoonanGenome.bed -b /home/ak2267/project/EnhancerClasses/'+Lift2_Dict[species]+'/me2/merge_*_overlap*_named.bed >> intermediate_'+species+'/'+species+'_jobNum'+str(i)+'_allIntersects.bed')
            outArray.append('bedtools intersect -wa -a intermediate_'+species+'/'+species+'_shuffled_jobNum'+str(i)+'_mapToNoonan_inNoonanGenome.bed -b /home/ak2267/project/EnhancerClasses/Limb/'+Lift2_Dict[species]+'/merge_*_overlap_named.bed >> intermediate_'+species+'/'+species+'_jobNum'+str(i)+'_allIntersects.bed')
        
            # count number with ANY mark
            outArray.append('allPeaks=$(cut -f 4 intermediate_'+species+'/'+species+'_jobNum'+str(i)+'_allIntersects.bed | sort -u | wc -l)')
            outArray.append('echo -e $allPeaks\"\\t\"$denominator >> countsWithPeak/'+species+'_all_all_countsWithPeak_jobNum'+str(i)+'.txt')
            
            # close for loop
            outArray.append('done')
            
        # write to CGI-centric job file 2
        outJobFile2.write(' ; '.join(outArray) +'\n')

# set number of jobs to write, and number of permutations per job
totalPermutations = 20000
numberOfJobs = 10
permutationsPerJob = totalPermutations / numberOfJobs

# write jobs for hg19 (much simpler)
for i in range(1,numberOfJobs+1):
    outArray = []
    species = 'hg19'

    # cd and source
    outArray.append('cd /gpfs/gibbs/pi/noonan/ak2267/singleSpecies/shuffling ; source /home/ak2267/.bashrc ; source /home/ak2267/.bash_profile')
    
    # calculate total oCGIs
    outArray.append('denominator=$(cat '+species+'_oCGIs_noFeatures.bed | wc -l)')
        
    # shuffle oCGIs
    outArray.append('for i in {1..'+str(permutationsPerJob)+'}; do bedtools shuffle -excl '+species+'_merge.bed -chrom -noOverlapping -i '+species+'_oCGIs_noFeatures.bed -g /home/ak2267/genomes/chrom.sizes/'+species+'.chrom.sizes > intermediate_'+species+'/'+species+'_shuffled_jobNum'+str(i)+'.bed')
    # note that I'm breaking here in the middle of a unix for loop
    
    # remove previously used file for collecting oCGIs with ANY peak to start fresh
    outArray.append('rm intermediate_'+species+'/'+species+'_jobNum'+str(i)+'_allIntersects.bed')
    
    # for each tissue x mark, intersect with peaks, and write numbers to countsWithPeak/${species}_${tissue}_${mark}_countsWithPeak_jobNum#.txt
    outArray.append('withPeak=$(bedtools intersect -wa -u -a intermediate_'+species+'/'+species+'_shuffled_jobNum'+str(i)+'.bed -b /home/ak2267/project/EnhancerClasses/'+species+'/ac/merge_*_overlap*_named.bed | wc -l)')
    outArray.append('echo -e $withPeak\"\\t\"$denominator >> countsWithPeak/'+species+'_devBrain_H3K27ac_countsWithPeak_jobNum'+str(i)+'.txt')
    outArray.append('withPeak=$(bedtools intersect -wa -u -a intermediate_'+species+'/'+species+'_shuffled_jobNum'+str(i)+'.bed -b /home/ak2267/project/EnhancerClasses/'+species+'/me2/merge_*_overlap*_named.bed | wc -l)')
    outArray.append('echo -e $withPeak\"\\t\"$denominator >> countsWithPeak/'+species+'_devBrain_H3K4me2_countsWithPeak_jobNum'+str(i)+'.txt')
    outArray.append('withPeak=$(bedtools intersect -wa -u -a intermediate_'+species+'/'+species+'_shuffled_jobNum'+str(i)+'.bed -b /home/ak2267/project/EnhancerClasses/Limb/'+species+'/merge_*_overlap_named.bed | wc -l)')
    outArray.append('echo -e $withPeak\"\\t\"$denominator >> countsWithPeak/'+species+'_devLimb_H3K27ac_countsWithPeak_jobNum'+str(i)+'.txt')
        
    # write to files for ANY mark calculations
    outArray.append('bedtools intersect -wa -a intermediate_'+species+'/'+species+'_shuffled_jobNum'+str(i)+'.bed -b /home/ak2267/project/EnhancerClasses/'+species+'/ac/merge_*_overlap*_named.bed >> intermediate_'+species+'/'+species+'_jobNum'+str(i)+'_allIntersects.bed')
    outArray.append('bedtools intersect -wa -a intermediate_'+species+'/'+species+'_shuffled_jobNum'+str(i)+'.bed -b /home/ak2267/project/EnhancerClasses/'+species+'/me2/merge_*_overlap*_named.bed >> intermediate_'+species+'/'+species+'_jobNum'+str(i)+'_allIntersects.bed')
    outArray.append('bedtools intersect -wa -a intermediate_'+species+'/'+species+'_shuffled_jobNum'+str(i)+'.bed -b /home/ak2267/project/EnhancerClasses/Limb/'+species+'/merge_*_overlap_named.bed >> intermediate_'+species+'/'+species+'_jobNum'+str(i)+'_allIntersects.bed')
    
    # count number with ANY mark
    outArray.append('allPeaks=$(cut -f 4 intermediate_'+species+'/'+species+'_jobNum'+str(i)+'_allIntersects.bed | sort -u | wc -l)')
    outArray.append('echo -e $allPeaks\"\\t\"$denominator >> countsWithPeak/'+species+'_all_all_countsWithPeak_jobNum'+str(i)+'.txt')
        
    # close for loop
    outArray.append('done')
            
    # write to CGI-centric job file 2
    outJobFile2.write(' ; '.join(outArray) +'\n')


##########################
###### PEAK-CENTRIC ######
##########################

# Job file 1: generate PEAKS in each species that reciprocally lift to human
# Job file 2:
# reshuffle onto same genome background
# lift to hg38 and back, only keeping ones that lift
# count sites with oCGI in each tissue x mark and write to separate output file for each tissue x mark
# ALSO lift sites to rheMac2 or mm9 as necessary in order to overlap with marks in those genomes

# Then return to 220912_CGI_singleSpecies_activity.sh to summarize using summarizeReshuffling_peakCentric.py

RollerSpecies = ['rheMac10', 'calJac4', 'mm39', 'rn7', 'susScr11', 'canFam6', 'felCat9', 'equCab3']
NoonanSpecies = ['hg19', 'rheMac2', 'mm9']

outJobFile1_peakCentric = open('230116_peakCentric_step1_get_peaks_jobFile.txt', 'wt')
outJobFile2_peakCentric = open('230116_peakCentric_step2_do_shuffling_jobFile.txt', 'wt')

# JOB FILE 1:
# generate PEAKS in each species that reciprocally lift to human and don't overlap human features
# the resulting files are in format ${species}_${tissue}_${mark}_noFeaturePeaks.bed

for species in RollerSpecies:

    filesToExclude = []

    for tissue in ['brain', 'liver', 'muscle', 'testis']:
    
        for mark in ['H3K4me3', 'H3K27ac', 'H3K4me1']:
        
            fileStem = species+'_'+tissue+'_'+mark
            
            # start collecting commands
            outArray = []
            outArray.append('cd /gpfs/gibbs/pi/noonan/ak2267/singleSpecies/shuffling ; source /home/ak2267/.bashrc ; source /home/ak2267/.bash_profile')
            
            # remove features from peak set            
            outArray.append('bedtools intersect -v -a /gpfs/gibbs/pi/noonan/ak2267/Roller/peaks/intersection/'+fileStem+'_intersection.bed -b '+species+'_merge.bed | awk \'{ print $1\"\\t\"$2\"\\t\"$3\"\\tPeak_\"NR }\' > '+fileStem+'_peaks_noFeaturesInRoller.bed')
            
            # lift to human and filter human annotations
            outArray.append('liftOver -minMatch=0.3 '+fileStem+'_peaks_noFeaturesInRoller.bed /home/ak2267/genomes/chain/'+species+'ToHg38.over.chain.gz '+fileStem+'_peaks_noFeaturesInRoller_LOhg38.bed unMapped')
            outArray.append('bedtools intersect -v -a '+fileStem+'_peaks_noFeaturesInRoller_LOhg38.bed -b hg38_merge.bed > '+fileStem+'_peaks_noFeaturesInRoller_LOhg38_noFeaturesInHg38.bed')
            outArray.append('liftOver -minMatch=0.3 '+fileStem+'_peaks_noFeaturesInRoller_LOhg38_noFeaturesInHg38.bed /home/ak2267/genomes/chain/hg38To'+species[0].upper()+species[1:]+'.over.chain.gz '+fileStem+'_peaks_noFeaturesInRoller_mapBack.bed unMapped')
            
            # filter to reciprocally lifting sites
            outArray.append('bedtools intersect -wao -a '+fileStem+'_peaks_noFeaturesInRoller_mapBack.bed -b '+fileStem+'_peaks_noFeaturesInRoller.bed > '+fileStem+'_peaks_noFeaturesInRoller_checkMapBack.txt')
            outArray.append('python /gpfs/gibbs/pi/noonan/ak2267/speciesPairs/consensusCGIs/Scripts/restrictToSitesThatMapBack.py '+fileStem+'_peaks_noFeaturesInRoller_checkMapBack.txt > '+fileStem+'_noFeaturePeaks.bed')
            # use this file downstream in step2 : fileStem+'_peaks_noFeatures.bed
            
            # write to output job file
            outJobFile1_peakCentric.write(' ; '.join(outArray)+'\n')
            

for species in NoonanSpecies:

    # write 3 lines for each species - brain ac, brain me2, limb ac
    FilePath_Dict = {}
    FilePath_Dict['brainAc'] = ['/home/ak2267/project/EnhancerClasses/', '/ac/merge_*_overlap*_named.bed']
    FilePath_Dict['brainMe2'] = ['/home/ak2267/project/EnhancerClasses/', '/me2/merge_*_overlap*_named.bed']
    FilePath_Dict['limbAc'] = ['/home/ak2267/project/EnhancerClasses/Limb/', '/merge_*_overlap_named.bed']
    
    if species == 'hg19':
    
        for tissueMark in FilePath_Dict:
            outJobFile1_peakCentric.write('cd /gpfs/gibbs/pi/noonan/ak2267/singleSpecies/shuffling ; source /home/ak2267/.bashrc ; source /home/ak2267/.bash_profile ; cat '+species.join(FilePath_Dict[tissueMark])+' > hg19_'+tissueMark+'_cat.bed ; bedtools intersect -v -a hg19_'+tissueMark+'_cat.bed -b hg19_merge.bed | sort -k1,1 -k2,2n - > hg19_'+tissueMark+'_sort.bed ; bedtools merge -i hg19_'+tissueMark+'_sort.bed | awk \'{ print $1\"\\t\"$2\"\\t\"$3\"\\tPeak_\"NR }\' > hg19_'+tissueMark+'_noFeaturePeaks.bed'+'\n')
        
    else:
        
        for tissueMark in FilePath_Dict:
            outArray = []
            outArray.append('cd /gpfs/gibbs/pi/noonan/ak2267/singleSpecies/shuffling ; source /home/ak2267/.bashrc ; source /home/ak2267/.bash_profile')
            
            fileStem = species + '_' + tissueMark
            
            # make merged list in rheMac2 / mm9
            outArray.append('cat '+species.join(FilePath_Dict[tissueMark])+' > '+fileStem+'_cat.bed ; bedtools intersect -v -a '+fileStem+'_cat.bed -b '+species+'_merge.bed | sort -k1,1 -k2,2n - > '+fileStem+'_sort.bed ; bedtools merge -i '+fileStem+'_sort.bed | awk \'{ print $1\"\\t\"$2\"\\t\"$3\"\\tPeak_\"NR }\' > '+fileStem+'_mergedPeaks.bed')
            
            # lift to human, restrict in human, lift back
            outArray.append('liftOver -minMatch=0.3 '+fileStem+'_mergedPeaks.bed /home/ak2267/genomes/chain/'+species+'ToHg19.over.chain.gz '+fileStem+'_mergedPeaks_LOhg19.bed unMapped')
            outArray.append('bedtools intersect -v -a '+fileStem+'_mergedPeaks_LOhg19.bed -b hg19_merge.bed > '+fileStem+'_mergedPeaks_LOhg19_noFeatureInHg19.bed')
            outArray.append('liftOver -minMatch=0.3 '+fileStem+'_mergedPeaks_LOhg19_noFeatureInHg19.bed /home/ak2267/genomes/chain/hg19To'+species[0].upper()+species[1:]+'.over.chain.gz '+fileStem+'_mergedPeaks_mapBack.bed unMapped')
            
            # restrict to reciprocally lifting sites and add name in col 4
            outArray.append('bedtools intersect -wao -a '+fileStem+'_mergedPeaks_mapBack.bed -b '+fileStem+'_mergedPeaks.bed > '+fileStem+'_mergedPeaks_checkMapBack.txt')
            outArray.append('python /gpfs/gibbs/pi/noonan/ak2267/speciesPairs/consensusCGIs/Scripts/restrictToSitesThatMapBack.py '+fileStem+'_mergedPeaks_checkMapBack.txt > '+fileStem+'_noFeaturePeaks.bed')
      
            # write to job file
            outJobFile1_peakCentric.write(' ; '.join(outArray)+'\n')
            
outJobFile1_peakCentric.close()


# JOB FILE 2: shuffling - these are the steps:
# 1) bedtools shuffle the PEAKS in each Roller species - exclude intervals in ${species}_merge.bed
# 2) lift these new intervals to hg38 and back, keeping only those that reciprocally lift AND do not have feature in hg38
# 3) for each tissue x mark, intersect this shuffle set restricted as described in step 2 with oCGIs
# 4) write number (numerator and denominator) to countsWithCGI/${species}_${tissue}_${mark}_countsWithCGI.txt

# same steps for Noonan genomes (rheMac2 / mm9 via hg19 lifts)
# 1) bedtools shuffle the PEAKS in each Noonan species - exclude intervals in ${species}_merge.bed
# 2) lift these new intervals to hg19 and back, keeping only those that reciprocally lift AND do not have feature in hg19
# 3) for each tissue x mark, intersect this shuffle set restricted as described in step 2 with oCGIs
# 4) write number (numerator and denominator) to countsWithCGI/${species}_${tissue}_${mark}_countsWithCGI.txt

# for hg19 peaks: simply shuffle within hg19

# For all species, do 20000 shuffles but split them into x separate jobs writing into y separate intermediate files (x * y = 20000)

# ROLLER SPECIES

for species in RollerSpecies:

    # set number of jobs to write, and number of permutations per job
    totalPermutations = 20000
    
    if species in ['susScr11', 'felCat9']:
        numberOfJobs = 80
        
    elif species in ['canFam6', 'equCab3']:
        numberOfJobs = 40
        
    else:
        numberOfJobs = 20
        
    permutationsPerJob = totalPermutations / numberOfJobs

    # loop through tissues
    for tissue in ['brain', 'liver', 'muscle', 'testis']:
        
        # loop through marks
        for mark in ['H3K4me3', 'H3K27ac', 'H3K4me1']:
        
            fileStem = species +'_'+ tissue +'_'+ mark
            
            # do this repeatedly for numberOfJobs times - use i to keep intermediate files distinct
            for i in range(1,numberOfJobs+1):
    
                outArray = []
    
                # cd and source
                outArray.append('cd /gpfs/gibbs/pi/noonan/ak2267/singleSpecies/shuffling ; source /home/ak2267/.bashrc ; source /home/ak2267/.bash_profile')
            
                # shuffle PEAKS
                outArray.append('for i in {1..'+str(permutationsPerJob)+'}; do bedtools shuffle -excl '+species+'_merge.bed -chrom -noOverlapping -i '+fileStem+'_noFeaturePeaks.bed -g /home/ak2267/genomes/chrom.sizes/'+species+'.chrom.sizes > peakIntermediate_'+species+'/'+fileStem+'_shuffled_jobNum'+str(i)+'.bed')
                # note that I'm breaking here in the middle of a unix for loop
        
                # lift to hg38 and back, with step in hg38 where I restrict on hg38 features
                outArray.append('liftOver -minMatch=0.3 peakIntermediate_'+species+'/'+fileStem+'_shuffled_jobNum'+str(i)+'.bed /home/ak2267/genomes/chain/'+species+'ToHg38.over.chain.gz peakIntermediate_'+species+'/'+fileStem+'_shuffled_jobNum'+str(i)+'_LOhg38.bed unMapped')
                outArray.append('bedtools intersect -v -a peakIntermediate_'+species+'/'+fileStem+'_shuffled_jobNum'+str(i)+'_LOhg38.bed -b hg38_merge.bed > peakIntermediate_'+species+'/'+fileStem+'_shuffled_jobNum'+str(i)+'_LOhg38_noFeatures.bed')
                outArray.append('liftOver -minMatch=0.3 peakIntermediate_'+species+'/'+fileStem+'_shuffled_jobNum'+str(i)+'_LOhg38_noFeatures.bed /home/ak2267/genomes/chain/hg38To'+species[0].upper()+species[1:]+'.over.chain.gz peakIntermediate_'+species+'/'+fileStem+'_shuffled_jobNum'+str(i)+'_LOhg38_mapBack.bed unMapped')
        
                # filter list in RollerSpecies genome to those that lift reciprocally
                outArray.append('bedtools intersect -wao -a peakIntermediate_'+species+'/'+fileStem+'_shuffled_jobNum'+str(i)+'_LOhg38_mapBack.bed -b peakIntermediate_'+species+'/'+fileStem+'_shuffled_jobNum'+str(i)+'.bed > peakIntermediate_'+species+'/'+fileStem+'_shuffled_jobNum'+str(i)+'_checkMapBack.txt')
                outArray.append('python /gpfs/gibbs/pi/noonan/ak2267/speciesPairs/consensusCGIs/Scripts/restrictToSitesThatMapBack.py peakIntermediate_'+species+'/'+fileStem+'_shuffled_jobNum'+str(i)+'_checkMapBack.txt > peakIntermediate_'+species+'/'+fileStem+'_shuffled_jobNum'+str(i)+'_mapBothWays.bed')
        
                # get denominator of sites that lifted reciprocally
                outArray.append('denominator=$(cat peakIntermediate_'+species+'/'+fileStem+'_shuffled_jobNum'+str(i)+'_mapBothWays.bed | wc -l)')
        
                # intersect with oCGIs, and write numbers to countsWithCGI/${species}_${tissue}_${mark}_countsWithCGI_jobNum#.txt
                outArray.append('withCGI=$(bedtools intersect -wa -u -a peakIntermediate_'+species+'/'+fileStem+'_shuffled_jobNum'+str(i)+'_mapBothWays.bed -b '+species+'_oCGIs_noFeatures.bed | wc -l)')
                outArray.append('echo -e $withCGI\"\\t\"$denominator >> countsWithCGI/'+fileStem+'_countsWithCGI_jobNum'+str(i)+'.txt')
            
                # close for loop
                outArray.append('done')
            
                # write to peak-centric job file 2
                outJobFile2_peakCentric.write(' ; '.join(outArray) +'\n')

# NOONAN SPECIES

tissueMarkArray = ['brainAc', 'brainMe2', 'limbAc']

for species in NoonanSpecies:

     # set number of jobs to write, and number of permutations per job
    totalPermutations = 20000
    
    if species in ['rheMac2', 'mm9']:
        numberOfJobs = 40
    else:
        numberOfJobs = 10
        
    permutationsPerJob = totalPermutations / numberOfJobs

    # loop through tissueMarks
    for tissueMark in tissueMarkArray:
        
        fileStem = species + '_' + tissueMark

        # separate job files (20000 / numberOfJobs)
        for i in range(1,numberOfJobs+1):
        
            outArray = []

            # cd and source
            outArray.append('cd /gpfs/gibbs/pi/noonan/ak2267/singleSpecies/shuffling ; source /home/ak2267/.bashrc ; source /home/ak2267/.bash_profile')
        
            # shuffle PEAKS
            outArray.append('for i in {1..'+str(permutationsPerJob)+'}; do bedtools shuffle -excl '+species+'_merge.bed -chrom -noOverlapping -i '+fileStem+'_noFeaturePeaks.bed -g /home/ak2267/genomes/chrom.sizes/'+species+'.chrom.sizes > peakIntermediate_'+species+'/'+fileStem+'_shuffled_jobNum'+str(i)+'.bed')
            # note that I'm breaking here in the middle of a unix for loop
            
            # if this is hg19, don't need to do any lifting steps - just count and close out the for loop
            if species == 'hg19':
            
                # calculate total PEAKS
                outArray.append('denominator=$(cat '+fileStem+'_noFeaturePeaks.bed | wc -l)')
                
                # intersect with oCGIs and write number to ${species}_${tissue}_${mark}_countsWithCGI.txt
                outArray.append('withCGI=$(bedtools intersect -wa -u -a peakIntermediate_'+species+'/'+fileStem+'_shuffled_jobNum'+str(i)+'.bed -b '+species+'_oCGIs_noFeatures.bed | wc -l)')
                outArray.append('echo -e $withCGI\"\\t\"$denominator >> countsWithCGI/'+fileStem+'_countsWithCGI_jobNum'+str(i)+'.txt')
                
                # close for loop
                outArray.append('done')
                
                # write to CGI-centric job file 2
                outJobFile2_peakCentric.write(' ; '.join(outArray) +'\n')
                
            # if this is rheMac2 or mm9, need to lift to hg19 and further restrict based on that
            if species == 'rheMac2' or species == 'mm9':
            
                # lift to hg19 and back, with step in hg19 where I restrict on hg19 features
                outArray.append('liftOver -minMatch=0.3 peakIntermediate_'+species+'/'+fileStem+'_shuffled_jobNum'+str(i)+'.bed /home/ak2267/genomes/chain/'+species+'ToHg19.over.chain.gz peakIntermediate_'+species+'/'+fileStem+'_shuffled_jobNum'+str(i)+'_LOhg19.bed unMapped')
                outArray.append('bedtools intersect -v -a peakIntermediate_'+species+'/'+fileStem+'_shuffled_jobNum'+str(i)+'_LOhg19.bed -b hg19_merge.bed > peakIntermediate_'+species+'/'+fileStem+'_shuffled_jobNum'+str(i)+'_LOhg19_noFeatures.bed')
                outArray.append('liftOver -minMatch=0.3 peakIntermediate_'+species+'/'+fileStem+'_shuffled_jobNum'+str(i)+'_LOhg19_noFeatures.bed /home/ak2267/genomes/chain/hg19To'+species[0].upper()+species[1:]+'.over.chain.gz peakIntermediate_'+species+'/'+fileStem+'_shuffled_jobNum'+str(i)+'_LOhg19_mapBack.bed unMapped')
        
                # filter list in NoonanSpecies genome to those that lift reciprocally
                outArray.append('bedtools intersect -wao -a peakIntermediate_'+species+'/'+fileStem+'_shuffled_jobNum'+str(i)+'_LOhg19_mapBack.bed -b peakIntermediate_'+species+'/'+fileStem+'_shuffled_jobNum'+str(i)+'.bed > peakIntermediate_'+species+'/'+fileStem+'_shuffled_jobNum'+str(i)+'_checkMapBack.txt')
                outArray.append('python /gpfs/gibbs/pi/noonan/ak2267/speciesPairs/consensusCGIs/Scripts/restrictToSitesThatMapBack.py peakIntermediate_'+species+'/'+fileStem+'_shuffled_jobNum'+str(i)+'_checkMapBack.txt > peakIntermediate_'+species+'/'+fileStem+'_shuffled_jobNum'+str(i)+'_mapBothWays.bed')

                # calculate total PEAKS
                outArray.append('denominator=$(cat peakIntermediate_'+species+'/'+fileStem+'_shuffled_jobNum'+str(i)+'_mapBothWays.bed | wc -l)')
                
                # intersect with oCGIs and write number to ${species}_${tissue}_${mark}_countsWithCGI.txt
                outArray.append('withCGI=$(bedtools intersect -wa -u -a peakIntermediate_'+species+'/'+fileStem+'_shuffled_jobNum'+str(i)+'.bed -b '+species+'_oCGIs_noFeatures.bed | wc -l)')
                outArray.append('echo -e $withCGI\"\\t\"$denominator >> countsWithCGI/'+fileStem+'_countsWithCGI_jobNum'+str(i)+'.txt')
                
                # close for loop
                outArray.append('done')
                    
                # write to CGI-centric job file 2
                outJobFile2_peakCentric.write(' ; '.join(outArray) +'\n')










