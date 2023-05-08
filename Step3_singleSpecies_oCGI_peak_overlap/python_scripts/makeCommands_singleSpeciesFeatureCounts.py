# 9/20/22
# Purpose: make commands file for running featureCounts on Roller peaks
# OR running bigWigAverageOverBed for Noonan peaks
# pulling commands from STEP6 in prepFiles from pairwise pipeline

# python makeCommands_singleSpeciesFeatureCounts.py /gpfs/gibbs/pi/noonan/ak2267/Roller/peaks/repsToUseRoller.txt > 220920_featureCounts_jobFile.txt

import sys

# ROLLER
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

for species in speciesArray:
    for tissue in Tissues:
        for mark in Marks:
        
    		#########
            # STEP 0: append info on where to work and source files
            outputString = []
            outputString.append('cd /gpfs/gibbs/pi/noonan/ak2267/singleSpecies/220919_activitySummary/featureCounts ; source /home/ak2267/.bashrc ; source /home/ak2267/.bash_profile')

            #STEP 6: run featureCounts to count reads in various peak intervals

            #turn into GTFs
            outputString.append('python /home/ak2267/Scripts/makeGTF_nameFromCol4.py /gpfs/gibbs/pi/noonan/ak2267/Roller/peaks/intersection/'+species+'_'+tissue+'_'+mark+'_intersection.bed '+species+'_'+tissue+'_'+mark+'_intersection.gtf')
            outputString.append('module load Subread/2.0.0-GCC-7.3.0-2.30')
            
            #run featureCounts
            fileList = []
            for replicate in Rep_Dict[species+'_'+tissue+'_'+mark]:
               outputString.append('featureCounts -t exon -g gene_id -a '+species+'_'+tissue+'_'+mark+'_intersection.gtf -o '+species+'_'+tissue+'_'+mark+'_'+replicate+'.counts /gpfs/gibbs/pi/noonan/ak2267/Roller/bam/'+species+'_'+tissue+'_'+mark+'_'+str(replicate)+'_bowtie2_filtered.bam')
               fileList.append(species+'_'+tissue+'_'+mark+'_'+replicate+'.counts')

            # generate file with RPM for each reconciledPeak (simpler than doing this in the final python script)
            outputString.append('python /gpfs/gibbs/pi/noonan/ak2267/speciesPairs/Scripts/countRPM_general.py '+str(len(fileList))+' '+' '.join(fileList)+' > '+species+'_'+tissue+'_'+mark+'.signal')
            
            print(' ; '.join(outputString))
            
# NOONAN

Tissues = ['brain','limb']
speciesArray = ['hg19','rheMac2','mm9']

Marks_Dict = {}
Marks_Dict['brain'] = ['ac','me2']
Marks_Dict['limb'] = ['ac']

# store timePoints for comparison in dictionary:
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
RepNumber_Dict['limb']['rheMac2'] = [0,3,1,0] # 3 replicates are all named rep1, but for timePoints e31, e32, e33
RepNumber_Dict['limb']['mm9'] =[3,2,2,2]


for species in speciesArray:
    for tissue in Tissues:
        for mark in Marks_Dict[tissue]:
            for timePointIndex in [0,1,2,3]:
            
            	#########
            	# STEP 0: append info on where to work and source files
            	outputString = []
            	outputString.append('cd /gpfs/gibbs/pi/noonan/ak2267/singleSpecies/220919_activitySummary/featureCounts ; source /home/ak2267/.bashrc ; source /home/ak2267/.bash_profile')

                pathToPeaks = TimePoints_Dict[tissue][species][timePointIndex].join(PeakPath_Dict[tissue][species][mark])
                # these are still arrays because they need both timePoint and rep filled in
                pathToBigWig = BigWigPath_Dict[tissue][species][mark]
                
                #STEP 6: run bigWigAverageOverBed to get average signal in peak intervals - THIS DIFFERS FROM ROLLER PIPELINE
                
                #run bigWigAverageOverBed
                fileList = []
                
                if timePointIndex == 1 and tissue == 'limb' and species == 'rheMac2':
                	outputString.append('bigWigAverageOverBed '+pathToBigWig[0]+'e31'+pathToBigWig[1]+'1'+pathToBigWig[2]+' '+pathToPeaks+' rheMac2_limb_ac_1_rep1.counts')
                	outputString.append('bigWigAverageOverBed '+pathToBigWig[0]+'e32'+pathToBigWig[1]+'1'+pathToBigWig[2]+' '+pathToPeaks+' rheMac2_limb_ac_1_rep2.counts')
                	outputString.append('bigWigAverageOverBed '+pathToBigWig[0]+'e33'+pathToBigWig[1]+'1'+pathToBigWig[2]+' '+pathToPeaks+' rheMac2_limb_ac_1_rep3.counts')
                	fileList.append('rheMac2_limb_ac_1_rep1.counts')
                	fileList.append('rheMac2_limb_ac_1_rep2.counts')
                	fileList.append('rheMac2_limb_ac_1_rep3.counts')
                elif RepNumber_Dict[tissue][species][timePointIndex] != 0:
                	for replicate in range(1,RepNumber_Dict[tissue][species][timePointIndex]+1):
                		# fix naming difference between mouse peaks and mouse bigWigs
                		if species == 'mm9' and tissue == 'brain' and (timePointIndex == 2 or timePointIndex == 3):
                			outputString.append('bigWigAverageOverBed '+pathToBigWig[0]+'e'+TimePoints_Dict[tissue][species][timePointIndex]+pathToBigWig[1]+str(replicate)+pathToBigWig[2]+' '+pathToPeaks+' '+species+'_'+tissue+'_'+mark+'_'+str(timePointIndex)+'_rep'+str(replicate)+'.counts')
                			fileList.append(species+'_'+tissue+'_'+mark+'_'+str(timePointIndex)+'_rep'+str(replicate)+'.counts')

                		else:
	                		outputString.append('bigWigAverageOverBed '+pathToBigWig[0]+TimePoints_Dict[tissue][species][timePointIndex]+pathToBigWig[1]+str(replicate)+pathToBigWig[2]+' '+pathToPeaks+' '+species+'_'+tissue+'_'+mark+'_'+str(timePointIndex)+'_rep'+str(replicate)+'.counts')
                			fileList.append(species+'_'+tissue+'_'+mark+'_'+str(timePointIndex)+'_rep'+str(replicate)+'.counts')
                
                # generate file with averageSignal for each reconciledPeak (simpler than doing this in the final python script)
                if RepNumber_Dict[tissue][species][timePointIndex] != 0:
                    outputString.append('python /gpfs/gibbs/pi/noonan/ak2267/speciesPairs/Scripts/countBigWigSignal_general.py '+str(len(fileList))+' '+' '.join(fileList)+' > '+species+'_'+tissue+'_'+mark+'_'+str(timePointIndex)+'.signal')
                
                    # print commands
                    print(' ; '.join(outputString))



