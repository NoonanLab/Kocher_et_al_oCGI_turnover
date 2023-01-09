# 8/10/22
# Purpose: make jobFile for doing analysis to make full spreadsheet for speciesPairs, CGIcentric
# Where each line in the spreadsheet is a CGI that maps between the species
# This file only Liver TF data - see other scripts for Roller histone mod data and Noonan histone mod data

import sys

speciesArray = ['rheMac10','mm39','rn7','canFam6']
Tissues = ['liver']
Marks = ['CEBPA','FOXA1','HNF4A','HNF6','CTCF']

SpeciesCombos = []

inReps = open(sys.argv[1],'rt')
Rep_Dict = {}
# Rep_Dict[species_tissue_mark] = [1,2] # or subset if a rep is unused
for line in inReps:
    splitLine = line.strip().split()
    Rep_Dict[splitLine[0]+'_'+splitLine[1]+'_'+splitLine[2]] = splitLine[3].split(',')

# store filePaths and fileEnds
workingDir = '/gpfs/gibbs/pi/noonan/ak2267/speciesPairs/CGIcentric'
pathToConsensusCGIs = '/gpfs/gibbs/pi/noonan/ak2267/speciesPairs/consensusCGIs/'
# then: {speciesA}_{speciesB}/{speciesA}_{speciesB}_speciesA_reconciledCGI.bed
#    or {speciesA}_{speciesB}/{speciesA}_{speciesB}_speciesB_reconciledCGI.bed
pathToAnnotations = '/home/ak2267/genomes/RefSeq/featureAnnotations/'
# then: ${species}_allFeatures.bed / _1_proteinCoding_promoters_merged.bed _2_proteinCoding_exons_merged.bed _3_lncRNA_merged.bed _4_ncRNA_merged.bed _5_pseudogene_merged.bed _6_unknown_merged.bed
pathToCTCF = '/gpfs/gibbs/pi/noonan/ak2267/LiverTF/peaks/CTCF/intersection/'
# then: ${species}_liver_CTCF_intersection.bed
pathToLiverTF = '/gpfs/gibbs/pi/noonan/ak2267/LiverTF/peaks/intersection/'
# then: ${species}_liver_${TF}_intersection.bed

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
            if mark == 'CTCF':
                pathToPeaks = pathToCTCF
                pathToBam = '/gpfs/gibbs/pi/noonan/ak2267/LiverTF/bam/CTCF/' # then format ${species}_liver_${TF}_${rep}_filtered.bam NOTE DIFFERENT THAN BELOW
            else:
                pathToPeaks = pathToLiverTF
                pathToBam = '/gpfs/gibbs/pi/noonan/ak2267/LiverTF/bam/' # then format ${species}_${TF}_${rep}_filtered.bam

            outputString = []

            #########
            # STEP 0: append info on where to work and source files
            outputString.append('cd '+workingDir+' ; source /home/ak2267/.bashrc ; source /home/ak2267/.bash_profile')
            outputString.append('mkdir LiverTF_'+speciesCombo+'_'+tissue+'_'+mark+' ; cd LiverTF_'+speciesCombo+'_'+tissue+'_'+mark)

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
            fileEndArray = ['_1_proteinCoding_promoters_merged.bed','_2_proteinCoding_exons_merged.bed','_3_lncRNA_merged.bed','_4_ncRNA_merged.bed','_5_pseudogene_merged.bed','_6_unknown_merged.bed']
            featuresArrayA = []
            featuresArrayB = []
            for fileEnd in fileEndArray:
                featuresArrayA.append(pathToAnnotations+speciesA+fileEnd)
                featuresArrayB.append(pathToAnnotations+speciesB+fileEnd)
            
            outputString.append('bedtools intersect -wa -u -a '+pathToPeaks+speciesA+'_'+tissue+'_'+mark+'_intersection.bed -b '+' -b '.join(featuresArrayA)+' | bedtools intersect -v -a '+reconciledA+' -b - > speciesA_intermediate_CGIsnoPromPeak.bed')
            outputString.append('bedtools intersect -wa -u -a '+pathToPeaks+speciesB+'_'+tissue+'_'+mark+'_intersection.bed -b '+' -b '.join(featuresArrayB)+' | bedtools intersect -v -a '+reconciledB+' -b - > speciesB_intermediate_CGIsnoPromPeak.bed')
            outputString.append('python /gpfs/gibbs/pi/noonan/ak2267/speciesPairs/consensusCGIs/Scripts/restrictToLO.py speciesA_intermediate_CGIsnoPromPeak.bed speciesB_intermediate_CGIsnoPromPeak.bed > speciesA_reconciledCGI_noFeaturePeaks.bed')
            outputString.append('python /gpfs/gibbs/pi/noonan/ak2267/speciesPairs/consensusCGIs/Scripts/restrictToLO.py speciesB_intermediate_CGIsnoPromPeak.bed speciesA_intermediate_CGIsnoPromPeak.bed > speciesB_reconciledCGI_noFeaturePeaks.bed')
            outputString.append('python /gpfs/gibbs/pi/noonan/ak2267/speciesPairs/consensusCGIs/Scripts/restrictToLO.py '+reconciledHuman+' speciesA_reconciledCGI_noFeaturePeaks.bed > human_reconciledCGI_noFeaturePeaks.bed')
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
                if mark == 'CTCF':
                    bamFileA = speciesA+'_'+tissue+'_'+mark+'_'+str(replicate)+'_filtered.bam'
                else:
                    bamFileA = speciesA+'_'+mark+'_'+str(replicate)+'_filtered.bam'

                outputString.append('featureCounts -t exon -g gene_id -a forFeatureCounts_A.gtf -o speciesA_'+mark+'_'+str(replicate)+'.counts '+pathToBam+bamFileA)
                A_fileList.append('speciesA_'+mark+'_'+str(replicate)+'.counts')

            for replicate in Rep_Dict[speciesB+'_'+tissue+'_'+mark]:
                if mark == 'CTCF':
                    bamFileB = speciesB+'_'+tissue+'_'+mark+'_'+str(replicate)+'_filtered.bam'
                else:
                    bamFileB = speciesB+'_'+mark+'_'+str(replicate)+'_filtered.bam'

                outputString.append('featureCounts -t exon -g gene_id -a forFeatureCounts_B.gtf -o speciesB_'+mark+'_'+str(replicate)+'.counts '+pathToBam+bamFileB)
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
            
            pythonInputFiles.append('repeatFilter')

            # STEP 8: python script to analyze the above intersection files & output final counts
            #echo $species $tissue $mark
            outputString.append('python /gpfs/gibbs/pi/noonan/ak2267/speciesPairs/Scripts/countFinalOverlaps_speciesPairs_CGIcentric.py '+' '.join(pythonInputFiles)+' >> '+workingDir+'/LiverTF_summaryFiles/'+speciesCombo+'_'+tissue+'_'+mark+'.txt')

            print('; '.join(outputString))
            
