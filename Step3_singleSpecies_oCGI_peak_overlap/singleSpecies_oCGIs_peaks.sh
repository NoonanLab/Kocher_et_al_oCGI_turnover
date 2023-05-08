# 9/12/22

# Purpose: table with these fields (CGI-centric):
# CGI_name CGI_length CpG_num 
# CGI_phastCons CGI_maxPhastCons CGI_totalPhastCons
# CGI_phastBias CGI_phastBiasLength CGI_phastBiasOverlap

# plus for Roller tissues:
# K4me3_overlap K27ac_overlap K4me1_overlap 
# Length_K4me3 Length_K27ac Length_K4me1 
# RPKM_K4me3 RPKM_K27ac RPKM_K4me1

# plus for Noonan tissues (note RPKM has different interpretation because it's from bigWig signal):
# K27ac_overlap K4me2_overlap 
# Length_K27ac Length_K4me2 
# RPKM_K27ac RPKM_K4me2

# Or the inverse, for each tissue (Peak - should also do this):
# Merged_peak_name Length_K4me3 Length_K27ac Length_K4me1 Length_union CGI CGI_length CpG_numCGI CpG_numPeak RPKM_K4me3 RPKM_K27ac RPKM_K4me1

# Then these can be merged across tissues in R

# work here
mkdir /gpfs/gibbs/pi/noonan/ak2267/singleSpecies/220919_activitySummary
cd /gpfs/gibbs/pi/noonan/ak2267/singleSpecies/220919_activitySummary

# CpG islands are here
/home/ak2267/genomes/CGI/UCSC_AL/${species}_CGIsAL.bed
# also need hg38_LOhg19_CGIsAL.bed, rheMac10_LOrheMac2_CGIsAL.bed, mm39_LOmm9_CGIsAL.bed

# feature annotations are here
/home/ak2267/genomes/RefSeq/featureAnnotations/${species}_allFeatures.bed

# PEAK FILES are here
# Roller ChIP peaks and signal:
/gpfs/gibbs/pi/noonan/ak2267/Roller/peaks/intersection/${species}_${tissue}_${mark}_intersection.bed
# Noonan BRAIN ChIP peaks and signal:
/home/ak2267/project/EnhancerClasses/hg19/ac/merge_${timePoint}_overlap_named.bed # timePoint = CS16/CS23/F2F/F2O
/home/ak2267/project/EnhancerClasses/hg19/me2/merge_${timePoint}_me2_overlap_named.bed # timePoint = CS16/CS23/F2F/F2O
/home/ak2267/project/EnhancerClasses/rheMac2/ac/merge_${timePoint}_overlap_rh_named.bed # timePoint = e55/e79F/e79O
/home/ak2267/project/EnhancerClasses/rheMac2/me2/merge_${timePoint}_overlap_me2_rh_named.bed # timePoint = e55/e79F/e79O
/home/ak2267/project/EnhancerClasses/mm9/ac/merge_${timePoint}_overlap_mm_named.bed # timePoint = e11/e14/17F/17O
/home/ak2267/project/EnhancerClasses/mm9/me2/merge_${timePoint}_overlap_me2_mm_named.bed # timePoint = e11/e14/17F/17O
# Noonan LIMB ChIP peaks and signal:
/home/ak2267/project/EnhancerClasses/Limb/hg19/merge_${timePoint}_overlap_named.bed # timePoint = E33/E41/E44/E47
/home/ak2267/project/EnhancerClasses/Limb/rheMac2/merge_${timePoint}_overlap_named.bed # timePoint = e31/e36
/home/ak2267/project/EnhancerClasses/Limb/mm9/merge_${timePoint}_overlap_named.bed # timePoint = e10.5/e11.5/e12.5/e13.5


################## CGI-centric
mkdir intersectPeaks
mkdir featureCounts
mkdir faCount
mkdir intersectPhastCons
mkdir intersectPhastBias

##### STEP 0: get intronic-intergenic CGIs

# restrict based on annotations in each species and add names in col 4
for species in rheMac10 calJac4 rn7 susScr11 canFam6 felCat9 equCab3
do
	echo $species
	bedtools intersect -v -a /home/ak2267/genomes/CGI/UCSC_AL/${species}_CGIsAL.bed \
	-b /home/ak2267/genomes/RefSeq/featureAnnotations/${species}_allFeatures.bed \
	| awk '{ print $1"\t"$2"\t"$3"\tCGI_"NR }' - \
	> ${species}_CGIsAL_noFeatures.bed
done

for species in hg38 mm39
do
	echo $species
	bedtools intersect -v -a /home/ak2267/genomes/CGI/UCSC_AL/${species}_CGIsAL.bed \
	-b /home/ak2267/genomes/RefSeq/featureAnnotations/${species}_allFeatures.bed \
	-b /home/ak2267/genomes/FANTOM/FANTOM_TSS_${species}.bed \
	-b /home/ak2267/genomes/blacklist/blacklist_${species}.bed \
	| awk '{ print $1"\t"$2"\t"$3"\tCGI_"NR }' - \
	> ${species}_CGIsAL_noFeatures.bed
done

# make files in hg19/rheMac2/mm9 with CGIs from hg38/rheMac10/mm39 that lift to these genome versions
# ${hg19/rheMac2/mm9}_CGIsAL_noFeatures.bed
# CGIs are named as in hg38/rheMac10/mm9 versions so they can be used together

liftOver -minMatch=0.8 hg38_CGIsAL_noFeatures.bed /home/ak2267/genomes/chain/hg38ToHg19.over.chain.gz \
hg38_CGIsAL_noFeatures_LOhg19.bed unMapped
liftOver -minMatch=0.8 hg38_CGIsAL_noFeatures_LOhg19.bed /home/ak2267/genomes/chain/hg19ToHg38.over.chain.gz \
hg38_CGIsAL_noFeatures_LOhg19_mapBack.bed unMapped
bedtools intersect -wao -a hg38_CGIsAL_noFeatures_LOhg19_mapBack.bed \
-b hg38_CGIsAL_noFeatures.bed > hg38_checkMapBack.txt
python /gpfs/gibbs/pi/noonan/ak2267/speciesPairs/consensusCGIs/Scripts/restrictToSitesThatMapBack.py \
hg38_checkMapBack.txt > hg38_CGIsAL_mapBothWays_inHg38.bed
python /gpfs/gibbs/pi/noonan/ak2267/speciesPairs/consensusCGIs/Scripts/restrictToLO.py \
hg38_CGIsAL_noFeatures_LOhg19.bed hg38_CGIsAL_mapBothWays_inHg38.bed \
> hg38_CGIsAL_noFeatures_LOhg19_uniqueMapping.bed
cut -f 1,2,3 /home/ak2267/genomes/CGI/UCSC_AL/hg19_CGIsAL.bed \
| bedtools intersect -wo -a - -b hg38_CGIsAL_noFeatures_LOhg19_uniqueMapping.bed \
| cut -f 1,2,3,7 > hg19_CGIsAL_noFeatures.bed # 20362 / 24487

liftOver -minMatch=0.8 rheMac10_CGIsAL_noFeatures.bed /home/ak2267/genomes/chain/rheMac10ToRheMac8.over.chain.gz \
rheMac10_CGIsAL_noFeatures_LOrheMac8.bed unMapped
liftOver -minMatch=0.8 rheMac10_CGIsAL_noFeatures_LOrheMac8.bed /home/ak2267/genomes/chain/rheMac8ToRheMac2.over.chain.gz \
rheMac10_CGIsAL_noFeatures_LOrheMac2.bed unMapped
liftOver -minMatch=0.8 rheMac10_CGIsAL_noFeatures_LOrheMac2.bed /home/ak2267/genomes/chain/rheMac2ToRheMac8.over.chain.gz \
rheMac10_CGIsAL_noFeatures_LOrheMac2_mapBack_rheMac8.bed unMapped
echo 'cd /gpfs/gibbs/pi/noonan/ak2267/singleSpecies/220919_activitySummary ; source /home/ak2267/.bashrc ; source /home/ak2267/.bash_profile ; liftOver -minMatch=0.8 rheMac10_CGIsAL_noFeatures_LOrheMac2_mapBack_rheMac8.bed /home/ak2267/genomes/chain/rheMac8ToRheMac10.over.chain.gz rheMac10_CGIsAL_noFeatures_LOrheMac2_mapBack.bed unMapped' > jobFile.txt
dsq --job-file jobFile.txt --mem-per-cpu 10G -c 1 --mail-type FAIL,END
sbatch dsq-jobFile-2022-09-19.sh # 21546946
bedtools intersect -wao -a rheMac10_CGIsAL_noFeatures_LOrheMac2_mapBack.bed \
-b rheMac10_CGIsAL_noFeatures.bed > rheMac10_checkMapBack.txt
python /gpfs/gibbs/pi/noonan/ak2267/speciesPairs/consensusCGIs/Scripts/restrictToSitesThatMapBack.py \
rheMac10_checkMapBack.txt > rheMac10_CGIsAL_mapBothWays_inRheMac10.bed
python /gpfs/gibbs/pi/noonan/ak2267/speciesPairs/consensusCGIs/Scripts/restrictToLO.py \
rheMac10_CGIsAL_noFeatures_LOrheMac2.bed rheMac10_CGIsAL_mapBothWays_inRheMac10.bed \
> rheMac10_CGIsAL_noFeatures_LOrheMac2_uniqueMapping.bed
cut -f 1,2,3 /home/ak2267/genomes/CGI/UCSC_AL/rheMac2_CGIsAL.bed \
| bedtools intersect -wo -a - -b rheMac10_CGIsAL_noFeatures_LOrheMac2_uniqueMapping.bed \
| cut -f 1,2,3,7 > rheMac2_CGIsAL_noFeatures.bed # 24480 / 38914

liftOver -minMatch=0.8 mm39_CGIsAL_noFeatures.bed /home/ak2267/genomes/chain/mm39ToMm10.over.chain.gz \
mm39_CGIsAL_noFeatures_LOmm10.bed unMapped
liftOver -minMatch=0.8 mm39_CGIsAL_noFeatures_LOmm10.bed /home/ak2267/genomes/chain/mm10ToMm9.over.chain.gz \
mm39_CGIsAL_noFeatures_LOmm9.bed unMapped
liftOver -minMatch=0.8 mm39_CGIsAL_noFeatures_LOmm9.bed /home/ak2267/genomes/chain/mm9ToMm10.over.chain.gz \
mm39_CGIsAL_noFeatures_LOmm9_mapBack_mm10.bed unMapped
liftOver -minMatch=0.8 mm39_CGIsAL_noFeatures_LOmm9_mapBack_mm10.bed /home/ak2267/genomes/chain/mm10ToMm39.over.chain.gz \
mm39_CGIsAL_noFeatures_LOmm9_mapBack.bed unMapped
bedtools intersect -wao -a mm39_CGIsAL_noFeatures_LOmm9_mapBack.bed \
-b mm39_CGIsAL_noFeatures.bed > mm39_checkMapBack.txt
python /gpfs/gibbs/pi/noonan/ak2267/speciesPairs/consensusCGIs/Scripts/restrictToSitesThatMapBack.py \
mm39_checkMapBack.txt > mm39_CGIsAL_mapBothWays_inMm39.bed
python /gpfs/gibbs/pi/noonan/ak2267/speciesPairs/consensusCGIs/Scripts/restrictToLO.py \
mm39_CGIsAL_noFeatures_LOmm9.bed mm39_CGIsAL_mapBothWays_inMm39.bed \
> mm39_CGIsAL_noFeatures_LOmm9_uniqueMapping.bed
cut -f 1,2,3 /home/ak2267/genomes/CGI/UCSC_AL/mm9_CGIsAL.bed \
| bedtools intersect -wo -a - -b mm39_CGIsAL_noFeatures_LOmm9_uniqueMapping.bed \
| cut -f 1,2,3,7 > mm9_CGIsAL_noFeatures.bed # 17252 / 18076
${species}_CGIsAL_noFeatures.bed
# NOTE THAT SOME GET SPLIT IN NOONAN GENOMES - NEED TO INCORPORATE MERGING IN PYTHON SCRIPT


##### STEP 1: intersect with ChIP-seq peaks

# Roller ChIP-seq (brain/liver/muscle/testis x H3K27ac/H3K4me3/H3K4me1)
for species in rheMac10 calJac4 mm39 rn7 susScr11 canFam6 felCat9 equCab3
do
for tissue in brain liver muscle testis
do
for mark in H3K27ac H3K4me3 H3K4me1
do
echo $species $tissue $mark
bedtools intersect -wao -a ${species}_CGIsAL_noFeatures.bed \
-b /gpfs/gibbs/pi/noonan/ak2267/Roller/peaks/intersection/${species}_${tissue}_${mark}_intersection.bed \
> intersectPeaks/${species}_Roller_${tissue}_${mark}_intersection.txt
done;done;done

# open for details on intersections with Noonan data
# Noonan brain ChIP-seq (4 tissue-timepoints x H3K27ac/H3K4me2)
for timePoint in CS16 CS23 F2F F2O
do
echo $timePoint
bedtools intersect -wao -a hg19_CGIsAL_noFeatures.bed \
-b /home/ak2267/project/EnhancerClasses/hg19/ac/merge_${timePoint}_overlap_named.bed \
> intersectPeaks/hg19_Noonan_brainAc_${timePoint}_intersection.txt
bedtools intersect -wao -a hg19_CGIsAL_noFeatures.bed \
-b /home/ak2267/project/EnhancerClasses/hg19/me2/merge_${timePoint}_me2_overlap_named.bed \
> intersectPeaks/hg19_Noonan_brainMe2_${timePoint}_intersection.txt
done

for timePoint in e55 e79F e79O
do
echo $timePoint
bedtools intersect -wao -a rheMac2_CGIsAL_noFeatures.bed \
-b /home/ak2267/project/EnhancerClasses/rheMac2/ac/merge_${timePoint}_overlap_rh_named.bed \
> intersectPeaks/rheMac2_Noonan_brainAc_${timePoint}_intersection.txt
bedtools intersect -wao -a rheMac2_CGIsAL_noFeatures.bed \
-b /home/ak2267/project/EnhancerClasses/rheMac2/me2/merge_${timePoint}_overlap_me2_rh_named.bed \
> intersectPeaks/rheMac2_Noonan_brainMe2_${timePoint}_intersection.txt
done

for timePoint in e11 e14 17F 17O
do
echo $timePoint
bedtools intersect -wao -a mm9_CGIsAL_noFeatures.bed \
-b /home/ak2267/project/EnhancerClasses/mm9/ac/merge_${timePoint}_overlap_mm_named.bed \
> intersectPeaks/mm9_Noonan_brainAc_${timePoint}_intersection.txt
bedtools intersect -wao -a mm9_CGIsAL_noFeatures.bed \
-b /home/ak2267/project/EnhancerClasses/mm9/me2/merge_${timePoint}_overlap_me2_mm_named.bed \
> intersectPeaks/mm9_Noonan_brainMe2_${timePoint}_intersection.txt
done

# Noonan limb ChIP-seq (4 tissue-timepoints x H3K27ac)
for timePoint in E33 E41 E44 E47
do
echo $timePoint
bedtools intersect -wao -a hg19_CGIsAL_noFeatures.bed \
-b /home/ak2267/project/EnhancerClasses/Limb/hg19/merge_${timePoint}_overlap_named.bed \
> intersectPeaks/hg19_Noonan_limbAc_${timePoint}_intersection.txt
done

for timePoint in e31 e36
do
echo $timePoint
bedtools intersect -wao -a rheMac2_CGIsAL_noFeatures.bed \
-b /home/ak2267/project/EnhancerClasses/Limb/rheMac2/merge_${timePoint}_overlap_named.bed \
> intersectPeaks/rheMac2_Noonan_limbAc_${timePoint}_intersection.txt
done

for timePoint in e10.5 e11.5 e12.5 e13.5
do
echo $timePoint
bedtools intersect -wao -a mm9_CGIsAL_noFeatures.bed \
-b /home/ak2267/project/EnhancerClasses/Limb/mm9/merge_${timePoint}_overlap_named.bed \
> intersectPeaks/mm9_Noonan_limbAc_${timePoint}_intersection.txt
done

# makes:
# ${species}_Noonan_brainAc_${timePoint}_intersection.txt
# ${species}_Noonan_brainMe2_${timePoint}_intersection.txt
# ${species}_Noonan_limbAc_${timePoint}_intersection.txt


##### STEP 2: run featureCounts to count reads in all peaks
cd /gpfs/gibbs/pi/noonan/ak2267/singleSpecies/220919_activitySummary

# requires python script to set up commands due to only using some replicates for Roller
cd /gpfs/gibbs/pi/noonan/ak2267/singleSpecies/220919_activitySummary/
python makeCommands_singleSpeciesFeatureCounts.py /gpfs/gibbs/pi/noonan/ak2267/Roller/peaks/repsToUseRoller.txt > 220920_featureCounts_jobFile.txt

dsq --job-file 220920_featureCounts_jobFile.txt --mem-per-cpu 5G -c 1 --mail-type FAIL,END
sbatch dsq-220920_featureCounts_jobFile-2022-09-22.sh # 17251376
sbatch dsq-220920_featureCounts_jobFile-2023-03-02.sh # 21546962

# makes these: 
/gpfs/gibbs/pi/noonan/ak2267/singleSpecies/220919_activitySummary/featureCounts/${species}_${tissue}_${mark}.signal				 # Roller
/gpfs/gibbs/pi/noonan/ak2267/singleSpecies/220919_activitySummary/featureCounts/${species}_${tissue}_${mark}_${timePoint}.signal # Noonan


##### STEP 3: count CpGs in CGIs and peaks using faCount
cd /gpfs/gibbs/pi/noonan/ak2267/singleSpecies/220919_activitySummary

# CGIs
for species in hg38 rheMac10 calJac4 mm39 rn7 susScr11 canFam6 felCat9 equCab3
do
echo 'cd /gpfs/gibbs/pi/noonan/ak2267/singleSpecies/220919_activitySummary; source /home/ak2267/.bashrc ; source /home/ak2267/.bash_profile ; bedtools getfasta -name -fi /gpfs/gibbs/pi/noonan/ak2267/genomes/'$species'.fa -bed '$species'_CGIsAL_noFeatures.bed > faCount/'$species'_CGIsAL.fa ; faCount faCount/'$species'_CGIsAL.fa > faCount/'$species'_CGIsAL.faCount'
done > 220920_CGIs_faCount_jobFile.txt

dsq --job-file 220920_CGIs_faCount_jobFile.txt --mem-per-cpu 5G -c 1 --mail-type FAIL,END
sbatch dsq-220920_CGIs_faCount_jobFile-2022-09-20.sh # 17141070
sbatch dsq-220920_CGIs_faCount_jobFile-2023-03-02.sh # 21547071
# re-run human (didn't have hg38.fa for first run) - only done for job 17141070, not necessary for 21547071
cd /gpfs/gibbs/pi/noonan/ak2267/singleSpecies/220919_activitySummary; source /home/ak2267/.bashrc ; source /home/ak2267/.bash_profile ; bedtools getfasta -name -fi /gpfs/gibbs/pi/noonan/ak2267/genomes/hg38.fa -bed hg38_CGIsAL_noFeatures.bed > faCount/hg38_CGIsAL.fa ; faCount faCount/hg38_CGIsAL.fa > faCount/hg38_CGIsAL.faCount

# ROLLER PEAKS
for species in rheMac10 calJac4 mm39 rn7 susScr11 canFam6 felCat9 equCab3
do
for tissue in brain liver muscle testis
do
for mark in H3K27ac H3K4me3 H3K4me1
do
echo 'cd /gpfs/gibbs/pi/noonan/ak2267/singleSpecies/220919_activitySummary; source /home/ak2267/.bashrc ; source /home/ak2267/.bash_profile ; bedtools getfasta -name -fi /gpfs/gibbs/pi/noonan/ak2267/genomes/'$species'.fa -bed /gpfs/gibbs/pi/noonan/ak2267/Roller/peaks/intersection/'$species'_'$tissue'_'$mark'_intersection.bed > faCount/'$species'_'$tissue'_'$mark'.fa ; faCount faCount/'$species'_'$tissue'_'$mark'.fa > faCount/'$species'_'$tissue'_'$mark'.faCount '
done;done;done > 220920_Roller_faCount_jobFile.txt

dsq --job-file 220920_Roller_faCount_jobFile.txt --mem-per-cpu 5G -c 1 --mail-type FAIL,END
sbatch dsq-220920_Roller_faCount_jobFile-2022-09-20.sh # 17140215
sbatch dsq-220920_Roller_faCount_jobFile-2023-03-02.sh # 21547113

# NOONAN - takes more lines because of timePoints being differently named across species

# Noonan brain
for timePoint in CS16 CS23 F2F F2O
do
echo 'cd /gpfs/gibbs/pi/noonan/ak2267/singleSpecies/220919_activitySummary; source /home/ak2267/.bashrc ; source /home/ak2267/.bash_profile ; bedtools getfasta -name -fi /gpfs/gibbs/pi/noonan/ak2267/genomes/hg19.fa -bed /home/ak2267/project/EnhancerClasses/hg19/ac/merge_'$timePoint'_overlap_named.bed > faCount/hg19_brain_ac_'$timePoint'.fa ; faCount faCount/hg19_brain_ac_'$timePoint'.fa > faCount/hg19_brain_ac_'$timePoint'.faCount'
echo 'cd /gpfs/gibbs/pi/noonan/ak2267/singleSpecies/220919_activitySummary; source /home/ak2267/.bashrc ; source /home/ak2267/.bash_profile ; bedtools getfasta -name -fi /gpfs/gibbs/pi/noonan/ak2267/genomes/hg19.fa -bed /home/ak2267/project/EnhancerClasses/hg19/me2/merge_'$timePoint'_me2_overlap_named.bed > faCount/hg19_brain_me2_'$timePoint'.fa ; faCount faCount/hg19_brain_me2_'$timePoint'.fa > faCount/hg19_brain_me2_'$timePoint'.faCount'
done > 220920_Noonan_faCount_jobFile.txt

for timePoint in e55 e79F e79O
do
echo 'cd /gpfs/gibbs/pi/noonan/ak2267/singleSpecies/220919_activitySummary; source /home/ak2267/.bashrc ; source /home/ak2267/.bash_profile ; bedtools getfasta -name -fi /gpfs/gibbs/pi/noonan/ak2267/genomes/rheMac2.fa -bed /home/ak2267/project/EnhancerClasses/rheMac2/ac/merge_'$timePoint'_overlap_rh_named.bed > faCount/rheMac2_brain_ac_'$timePoint'.fa ; faCount faCount/rheMac2_brain_ac_'$timePoint'.fa > faCount/rheMac2_brain_ac_'$timePoint'.faCount'
echo 'cd /gpfs/gibbs/pi/noonan/ak2267/singleSpecies/220919_activitySummary; source /home/ak2267/.bashrc ; source /home/ak2267/.bash_profile ; bedtools getfasta -name -fi /gpfs/gibbs/pi/noonan/ak2267/genomes/rheMac2.fa -bed /home/ak2267/project/EnhancerClasses/rheMac2/me2/merge_'$timePoint'_overlap_me2_rh_named.bed > faCount/rheMac2_brain_me2_'$timePoint'.fa ; faCount faCount/rheMac2_brain_me2_'$timePoint'.fa > faCount/rheMac2_brain_me2_'$timePoint'.faCount'
done >> 220920_Noonan_faCount_jobFile.txt

for timePoint in e11 e14 17F 17O
do
echo 'cd /gpfs/gibbs/pi/noonan/ak2267/singleSpecies/220919_activitySummary; source /home/ak2267/.bashrc ; source /home/ak2267/.bash_profile ; bedtools getfasta -name -fi /gpfs/gibbs/pi/noonan/ak2267/genomes/mm9.fa -bed /home/ak2267/project/EnhancerClasses/mm9/ac/merge_'$timePoint'_overlap_mm_named.bed > faCount/mm9_brain_ac_'$timePoint'.fa ; faCount faCount/mm9_brain_ac_'$timePoint'.fa > faCount/mm9_brain_ac_'$timePoint'.faCount'
echo 'cd /gpfs/gibbs/pi/noonan/ak2267/singleSpecies/220919_activitySummary; source /home/ak2267/.bashrc ; source /home/ak2267/.bash_profile ; bedtools getfasta -name -fi /gpfs/gibbs/pi/noonan/ak2267/genomes/mm9.fa -bed /home/ak2267/project/EnhancerClasses/mm9/me2/merge_'$timePoint'_overlap_me2_mm_named.bed > faCount/mm9_brain_me2_'$timePoint'.fa ; faCount faCount/mm9_brain_me2_'$timePoint'.fa > faCount/mm9_brain_me2_'$timePoint'.faCount'
done >> 220920_Noonan_faCount_jobFile.txt

# Noonan limb
for timePoint in E33 E41 E44 E47
do
echo 'cd /gpfs/gibbs/pi/noonan/ak2267/singleSpecies/220919_activitySummary; source /home/ak2267/.bashrc ; source /home/ak2267/.bash_profile ; bedtools getfasta -name -fi /gpfs/gibbs/pi/noonan/ak2267/genomes/hg19.fa -bed /home/ak2267/project/EnhancerClasses/Limb/hg19/merge_'$timePoint'_overlap_named.bed > faCount/hg19_limb_ac_'$timePoint'.fa ; faCount faCount/hg19_limb_ac_'$timePoint'.fa > faCount/hg19_limb_ac_'$timePoint'.faCount'
done >> 220920_Noonan_faCount_jobFile.txt

for timePoint in e31 e36
do
echo 'cd /gpfs/gibbs/pi/noonan/ak2267/singleSpecies/220919_activitySummary; source /home/ak2267/.bashrc ; source /home/ak2267/.bash_profile ; bedtools getfasta -name -fi /gpfs/gibbs/pi/noonan/ak2267/genomes/rheMac2.fa -bed /home/ak2267/project/EnhancerClasses/Limb/rheMac2/merge_'$timePoint'_overlap_named.bed > faCount/rheMac2_limb_ac_'$timePoint'.fa ; faCount faCount/rheMac2_limb_ac_'$timePoint'.fa > faCount/rheMac2_limb_ac_'$timePoint'.faCount'
done >> 220920_Noonan_faCount_jobFile.txt

for timePoint in e10.5 e11.5 e12.5 e13.5
do
echo 'cd /gpfs/gibbs/pi/noonan/ak2267/singleSpecies/220919_activitySummary; source /home/ak2267/.bashrc ; source /home/ak2267/.bash_profile ; bedtools getfasta -name -fi /gpfs/gibbs/pi/noonan/ak2267/genomes/mm9.fa -bed /home/ak2267/project/EnhancerClasses/Limb/mm9/merge_'$timePoint'_overlap_named.bed > faCount/mm9_limb_ac_'$timePoint'.fa ; faCount faCount/mm9_limb_ac_'$timePoint'.fa > faCount/mm9_limb_ac_'$timePoint'.faCount'
done >> 220920_Noonan_faCount_jobFile.txt

dsq --job-file 220920_Noonan_faCount_jobFile.txt --mem-per-cpu 5G -c 1 --mail-type FAIL,END
sbatch dsq-220920_Noonan_faCount_jobFile-2022-09-20.sh # 17140921
sbatch dsq-220920_Noonan_faCount_jobFile-2023-03-02.sh # 21547355

# makes these:
/gpfs/gibbs/pi/noonan/ak2267/singleSpecies/220919_activitySummary/faCount/${species}_CGIsAL.faCount			    			# CGIs
/gpfs/gibbs/pi/noonan/ak2267/singleSpecies/220919_activitySummary/faCount/${species}_${tissue}_${mark}.faCount			    # Roller
/gpfs/gibbs/pi/noonan/ak2267/singleSpecies/220919_activitySummary/faCount/${species}_${tissue}_${mark}_${timePoint}.faCount # Noonan



##### STEP 4: intersect phastCons with CGIs and peaks - tricky because phastCons is in hg38
# will need a column for "Lifts to hg38" which ALSO incorporates not being annotated as a feature in hg38

# lift all species to hg38 - will just use CGIs from Roller species and integrate Noonan species into that spreadsheet
# this takes a long time to run, should have run as jobs
cd /gpfs/gibbs/pi/noonan/ak2267/singleSpecies/220919_activitySummary/intersectPhastCons

for species in rheMac10 calJac4 mm39 rn7 susScr11 canFam6 felCat9 equCab3
do
echo $species
# to hg38
liftOver -minMatch=0.3 ../${species}_CGIsAL_noFeatures.bed /home/ak2267/genomes/chain/${species}ToHg38.over.chain.gz ${species}_CGIsAL_noFeatures_LOhg38.bed unMapped
# back to species
liftOver -minMatch=0.3 ${species}_CGIsAL_noFeatures_LOhg38.bed /home/ak2267/genomes/chain/hg38To${species^}.over.chain.gz ${species}_CGIsAL_noFeatures_LOhg38_mapBack.bed unMapped
# intersect & run restriction python script
bedtools intersect -wao -a ${species}_CGIsAL_noFeatures_LOhg38_mapBack.bed -b ../${species}_CGIsAL_noFeatures.bed > ${species}_checkMapBack.txt
python /gpfs/gibbs/pi/noonan/ak2267/speciesPairs/consensusCGIs/Scripts/restrictToSitesThatMapBack.py ${species}_checkMapBack.txt > ${species}_CGIsAL_mapBothWays_in${species}.bed
python /gpfs/gibbs/pi/noonan/ak2267/speciesPairs/consensusCGIs/Scripts/restrictToLO.py ${species}_CGIsAL_noFeatures_LOhg38.bed ${species}_CGIsAL_mapBothWays_in${species}.bed > ${species}_CGIsAL_noFeatures_LOhg38_uniqueMapping.bed
# remove sites with features in hg38
bedtools intersect -v -a ${species}_CGIsAL_noFeatures_LOhg38_uniqueMapping.bed -b /home/ak2267/genomes/RefSeq/featureAnnotations/hg38_allFeatures.bed -b /home/ak2267/genomes/blacklist/blacklist_hg38.bed -b /home/ak2267/genomes/FANTOM/FANTOM_TSS_hg38.bed > ${species}_CGIsAL_noFeatures_LOhg38_uniqueMapping_noHg38features.bed
# intersect with phastCons in hg38
bedtools intersect -wao -a ${species}_CGIsAL_noFeatures_LOhg38_uniqueMapping_noHg38features.bed -b /home/ak2267/genomes/phastCons/phastConsElements100way_hg38.bed > ${species}_CGIsAL_intersectPhastConsInHg38.txt
done
species=hg38
bedtools intersect -wao -a ../${species}_CGIsAL_noFeatures.bed -b /home/ak2267/genomes/phastCons/phastConsElements100way_hg38.bed > ${species}_CGIsAL_intersectPhastConsInHg38.txt

# makes:
/gpfs/gibbs/pi/noonan/ak2267/singleSpecies/220919_activitySummary/intersectPhastCons/${species}_CGIsAL_intersectPhastConsInHg38.txt
 
# removing sites annotated as features in hg38 removes a good chunk of sites (~20%)


##### STEP 5: intersect phastBias with CGIs and peaks - tricky because phastBias is in hg38
# use lifted files from above with features removed:
/gpfs/gibbs/pi/noonan/ak2267/singleSpecies/220919_activitySummary/intersectPhastCons/${species}_CGIsAL_noFeatures_LOhg38_uniqueMapping_noHg38features.bed

# phastBias files:
/gpfs/gibbs/pi/noonan/ak2267/speciesPairs/SequenceAnalysis/gBGC_analysis/phastBias_bed/${species}_allChr_phastBias_merged.bed

cd /gpfs/gibbs/pi/noonan/ak2267/singleSpecies/220919_activitySummary/intersectPhastBias

# intersect with phastBias file for that species (which are in hg38)
# have to alternate commands between consecutive species assignments (upper) and bedtools intersect (lower) to cover all species
species1=rheMac10; species2=rheMac8
species1=calJac4; species2=calJac3
species1=mm39; species2=mm10
species1=rn7; species2=rn6
species1=susScr11; species2=susScr11
species1=canFam6; species2=canFam3
species1=felCat9; species2=felCat8
species1=equCab3; species2=HLequCab3
bedtools intersect -wao -a /gpfs/gibbs/pi/noonan/ak2267/singleSpecies/220919_activitySummary/intersectPhastCons/${species1}_CGIsAL_noFeatures_LOhg38_uniqueMapping_noHg38features.bed -b /gpfs/gibbs/pi/noonan/ak2267/speciesPairs/SequenceAnalysis/gBGC_analysis/phastBias_bed/${species2}_allChr_phastBias_merged.bed > ${species1}_CGIsAL_intersectPhastBiasInHg38.txt

species=hg38
bedtools intersect -wao -a ../${species}_CGIsAL_noFeatures.bed -b /gpfs/gibbs/pi/noonan/ak2267/speciesPairs/SequenceAnalysis/gBGC_analysis/phastBias_bed/${species}_allChr_phastBias_merged.bed > ${species}_CGIsAL_intersectPhastBiasInHg38.txt

# makes:
/gpfs/gibbs/pi/noonan/ak2267/singleSpecies/220919_activitySummary/intersectPhastBias/${species}_CGIsAL_intersectPhastBiasInHg38.txt

##### STEP 6: integrate all info with python script
cd /gpfs/gibbs/pi/noonan/ak2267/singleSpecies/220919_activitySummary
mkdir /gpfs/gibbs/pi/noonan/ak2267/singleSpecies/220919_activitySummary/summaryTables_singleSpeciesCGIcentric

for species in hg38 rheMac10 calJac4 mm39 rn7 susScr11 canFam6 felCat9 equCab3
do
	echo $species
	python makeSingleSpeciesSummaryTables_CGIcentric.py ${species}
done

tar -zcvf summaryTables_singleSpeciesCGIcentric.gz summaryTables_singleSpeciesCGIcentric/
cd /Users/acadiak/Desktop/CGI/singleSpecies
scp ak2267@ruddle.hpc.yale.edu:/gpfs/gibbs/pi/noonan/ak2267/singleSpecies/220919_activitySummary/summaryTables_singleSpeciesCGIcentric.gz .



################# Peak-centric
# 9/29/22

# work here
mkdir /gpfs/gibbs/pi/noonan/ak2267/singleSpecies/220929_activitySummary_peakCentric
cd /gpfs/gibbs/pi/noonan/ak2267/singleSpecies/220929_activitySummary_peakCentric

mkdir featureAnnotations
mkdir liftToHuman
mkdir intersectCGIs
mkdir featureCounts
mkdir faCount
mkdir intersectPhastCons
mkdir intersectAge
mkdir intersectRepeats
mkdir summaryTables_singleSpeciesPeakCentric

# CpG islands are here
/home/ak2267/genomes/CGI/UCSC_AL/${species}_CGIsAL.bed
# if each tissue x mark x timePoint is separate, want to use CGIs called in each genome version 
# for hg19/rheMac2/mm9 instead of lifting

# feature annotations are here
/home/ak2267/genomes/RefSeq/featureAnnotations/${species}_allFeatures.bed
/home/ak2267/genomes/blacklist/blacklist_${species}.bed		# human and mouse only
/home/ak2267/genomes/FANTOM/FANTOM_TSS_${species}.bed		# human and mouse only

# PEAK FILES are here
# Roller ChIP peaks and signal:
/gpfs/gibbs/pi/noonan/ak2267/Roller/peaks/intersection/${species}_${tissue}_${mark}_intersection.bed
# Noonan BRAIN ChIP peaks and signal:
/home/ak2267/project/EnhancerClasses/hg19/ac/merge_${timePoint}_overlap_named.bed # timePoint = CS16/CS23/F2F/F2O
/home/ak2267/project/EnhancerClasses/hg19/me2/merge_${timePoint}_me2_overlap_named.bed # timePoint = CS16/CS23/F2F/F2O
/home/ak2267/project/EnhancerClasses/rheMac2/ac/merge_${timePoint}_overlap_rh_named.bed # timePoint = e55/e79F/e79O
/home/ak2267/project/EnhancerClasses/rheMac2/me2/merge_${timePoint}_overlap_me2_rh_named.bed # timePoint = e55/e79F/e79O
/home/ak2267/project/EnhancerClasses/mm9/ac/merge_${timePoint}_overlap_mm_named.bed # timePoint = e11/e14/17F/17O
/home/ak2267/project/EnhancerClasses/mm9/me2/merge_${timePoint}_overlap_me2_mm_named.bed # timePoint = e11/e14/17F/17O
# Noonan LIMB ChIP peaks and signal:
/home/ak2267/project/EnhancerClasses/Limb/hg19/merge_${timePoint}_overlap_named.bed # timePoint = E33/E41/E44/E47
/home/ak2267/project/EnhancerClasses/Limb/rheMac2/merge_${timePoint}_overlap_named.bed # timePoint = e31/e36
/home/ak2267/project/EnhancerClasses/Limb/mm9/merge_${timePoint}_overlap_named.bed # timePoint = e10.5/e11.5/e12.5/e13.5


##### STEP 0: collect ChIP-seq peaks in each species
# see file paths above
# Roller ChIP-seq (brain/liver/muscle/testis x H3K27ac/H3K4me3/H3K4me1)
# Noonan brain ChIP-seq (4 tissue-timepoints x H3K27ac/H3K4me2)
# Noonan limb ChIP-seq (4 tissue-timepoints x H3K27ac)

##### STEP 1: restrict to intronic-intergenic in each species
# then lift to human (with lift back to check specificity) and restrict to intronic-intergenic in human
# finally intersect with CGIs in original species

# do this and next step with a python script to make job files
cd /gpfs/gibbs/pi/noonan/ak2267/singleSpecies/220929_activitySummary_peakCentric/ 
python prepFiles_singleSpecies_peakCentric.py > 221017_singleSpecies_peakCentric_jobFile.txt

dsq --job-file 221017_singleSpecies_peakCentric_jobFile.txt --mem-per-cpu 5G -c 1 --mail-type FAIL,END
sbatch dsq-221017_singleSpecies_peakCentric_jobFile-2022-10-17.sh # 17881095
sbatch dsq-221017_singleSpecies_peakCentric_jobFile-2023-03-05.sh # 21635951

# gzip and download for use in R
tar -zcvf summaryFiles_singleSpeciesPeakCentric.gz summaryFiles_singleSpeciesPeakCentric/
cd /Users/acadiak/Desktop/CGI/singleSpecies/peakCentric
scp ak2267@ruddle.hpc.yale.edu:/gpfs/gibbs/pi/noonan/ak2267/singleSpecies/220929_activitySummary_peakCentric/summaryFiles_singleSpeciesPeakCentric.gz .





##########################################################################################
# Run reshuffling/permutation test for enrichment over intronic/intergenic background
# 1/11/23

mkdir /gpfs/gibbs/pi/noonan/ak2267/singleSpecies/shuffling
cd /gpfs/gibbs/pi/noonan/ak2267/singleSpecies/shuffling

##### General pipeline:
# Shuffle oCGIs on intronic-intergenic sequence within the species
# Lift shuffled oCGIs to human (with step to lift back, to ensure reciprocal mapping)
# Count the shuffled percentage with a peak = # that lift that overlap a peak / # that lift

##### Make regions to exclude from the reshuffling background

## Merge with features in each species to make bed file for exclusion
## Also remove rmsk sites from all EXCEPT hg38 which is only referenced after lifts 

# ROLLER GENOMES
for species in rheMac10 calJac4 rn7 susScr11 canFam6 felCat9 equCab3 # all except mm39
do
	echo $species
	cat /home/ak2267/genomes/RefSeq/featureAnnotations/${species}_allFeatures.bed /home/ak2267/genomes/rmsk/rmsk_${species}_merged.bed | cut -f 1,2,3 > ${species}_cat.bed
	sort -k1,1 -k2,2n ${species}_cat.bed > ${species}_sort.bed
	bedtools merge -i ${species}_sort.bed > ${species}_merge.bed
done

# add FANTOM and blacklist for mouse
species=mm39
cat /home/ak2267/genomes/RefSeq/featureAnnotations/${species}_allFeatures.bed /home/ak2267/genomes/FANTOM/FANTOM_TSS_mm39.bed /home/ak2267/genomes/blacklist/blacklist_mm39.bed /home/ak2267/genomes/rmsk/rmsk_${species}_merged.bed | cut -f 1,2,3 > ${species}_cat.bed
sort -k1,1 -k2,2n ${species}_cat.bed > ${species}_sort.bed
bedtools merge -i ${species}_sort.bed > ${species}_merge.bed

# NOONAN GENOMES
species=hg19
cat /home/ak2267/genomes/RefSeq/featureAnnotations/hg19_allFeatures.bed /home/ak2267/genomes/FANTOM/FANTOM_TSS_hg19.bed /home/ak2267/genomes/blacklist/blacklist_hg19.bed /home/ak2267/genomes/rmsk/rmsk_${species}_merged.bed | cut -f 1,2,3 > ${species}_cat.bed
sort -k1,1 -k2,2n ${species}_cat.bed > ${species}_sort.bed
bedtools merge -i ${species}_sort.bed > ${species}_merge.bed
species=rheMac2
cat /home/ak2267/genomes/RefSeq/featureAnnotations/${species}_allFeatures.bed /home/ak2267/genomes/rmsk/rmsk_${species}_merged.bed | cut -f 1,2,3 > ${species}_cat.bed
sort -k1,1 -k2,2n ${species}_cat.bed > ${species}_sort.bed
bedtools merge -i ${species}_sort.bed > ${species}_merge.bed
species=mm9
cat /home/ak2267/genomes/RefSeq/featureAnnotations/${species}_allFeatures.bed /home/ak2267/genomes/FANTOM/FANTOM_TSS_mm9.bed /home/ak2267/genomes/blacklist/blacklist_mm9.bed /home/ak2267/genomes/rmsk/rmsk_${species}_merged.bed | cut -f 1,2,3 > ${species}_cat.bed
sort -k1,1 -k2,2n ${species}_cat.bed > ${species}_sort.bed
bedtools merge -i ${species}_sort.bed > ${species}_merge.bed

rm *_cat.bed
rm *_sort.bed

# DO hg38 FOR REMOVING FEATURES AFTER LIFTS TO HUMAN - do NOT include repeatMasker track
species=hg38
cat /home/ak2267/genomes/RefSeq/featureAnnotations/hg38_allFeatures.bed /home/ak2267/genomes/FANTOM/FANTOM_TSS_hg38.bed /home/ak2267/genomes/blacklist/blacklist_hg38.bed | cut -f 1,2,3 > ${species}_cat.bed
sort -k1,1 -k2,2n ${species}_cat.bed > ${species}_sort.bed
bedtools merge -i ${species}_sort.bed > ${species}_merge.bed


####################################################
##### Run permutation test on CGI-centric data #####
####################################################

cd /gpfs/gibbs/pi/noonan/ak2267/singleSpecies/shuffling

# CGI file path: /home/ak2267/genomes/CGI/UCSC_AL/${species}_CGIsAL.bed

rm *_oCGIcounts.txt
rm *_countsWithPeak.txt
rm shuffleSummary.txt

# run prepFiles_singleSpecies_shuffling.py to make job files:
# 230112_step1_get_oCGIs_jobFile.txt
# 230112_step2_do_shuffling_jobFile.txt
python prepFiles_singleSpecies_shuffling.py

# run CGI-centric job file 1 (prepare bed files in Roller & Noonan genomes)
dsq --job-file 230112_step1_get_oCGIs_jobFile.txt --mem-per-cpu 5G -c 1 --mail-type FAIL,END
sbatch dsq-230112_step1_get_oCGIs_jobFile-2023-01-12.sh # 19899722
sbatch dsq-230112_step1_get_oCGIs_jobFile-2023-03-02.sh # 21576773

# collect OBSERVED VALUES into observed_counts.txt
for species in rheMac10 calJac4 mm39 rn7 susScr11 canFam6 felCat9 equCab3
do
    denominator=$(cat ${species}_oCGIs_noFeatures.bed | wc -l)
    for tissue in brain liver muscle testis
    do
        for mark in H3K4me3 H3K27ac H3K4me1
        do
            withPeak=$(bedtools intersect -wa -u -a ${species}_oCGIs_noFeatures.bed -b /gpfs/gibbs/pi/noonan/ak2267/Roller/peaks/intersection/${species}_${tissue}_${mark}_intersection.bed | wc -l)
            echo -e $species"\t"$tissue"\t"$mark"\t"$withPeak"\t"$denominator
        done; done; done >> observedCounts.txt

for species in hg19 rheMac2 mm9
do
    denominator=$(cat ${species}_oCGIs_noFeatures.bed | wc -l)
    
    tissue=devBrain
    for mark in ac me2
    do
        withPeak=$(bedtools intersect -wa -u -a ${species}_oCGIs_noFeatures.bed -b /home/ak2267/project/EnhancerClasses/${species}/${mark}/merge_*_overlap*_named.bed | wc -l)
        echo -e $species"\t"$tissue"\t"$mark"\t"$withPeak"\t"$denominator >> observedCounts.txt
    done
        
    tissue=devLimb
    mark=ac
    withPeak=$(bedtools intersect -wa -u -a ${species}_oCGIs_noFeatures.bed -b /home/ak2267/project/EnhancerClasses/Limb/${species}/merge_*_overlap_named.bed | wc -l)
    echo -e $species"\t"$tissue"\t"$mark"\t"$withPeak"\t"$denominator
done >> observedCounts.txt

# collect OBSERVED VALUES for ANY TISSUE AND MARK
for species in rheMac10 calJac4 mm39 rn7 susScr11 canFam6 felCat9 equCab3
do
    for tissue in brain liver muscle testis
    do
        for mark in H3K4me3 H3K27ac H3K4me1
        do
            bedtools intersect -wa -u -a ${species}_oCGIs_noFeatures.bed -b /gpfs/gibbs/pi/noonan/ak2267/Roller/peaks/intersection/${species}_${tissue}_${mark}_intersection.bed >> ${species}_allIntersects.bed
        done; done; done 

for species in calJac4 rn7 susScr11 canFam6 felCat9 equCab3
do
    denominator=$(cat ${species}_oCGIs_noFeatures.bed | wc -l)
    withPeak=$(cut -f 4 ${species}_allIntersects.bed | sort -u | wc -l)
    echo -e $species"\t"all"\t"all"\t"$withPeak"\t"$denominator
done >> observedCounts.txt

species=hg19
denominator=$(cat ${species}_oCGIs_noFeatures.bed | wc -l)
for mark in ac me2
do
    bedtools intersect -wa -u -a ${species}_oCGIs_noFeatures.bed -b /home/ak2267/project/EnhancerClasses/${species}/${mark}/merge_*_overlap*_named.bed >> ${species}_allIntersects.bed
done
mark=ac
bedtools intersect -wa -u -a ${species}_oCGIs_noFeatures.bed -b /home/ak2267/project/EnhancerClasses/Limb/${species}/merge_*_overlap_named.bed >> ${species}_allIntersects.bed
withPeak=$(cut -f 4 ${species}_allIntersects.bed | sort -u | wc -l)
echo -e $species"\t"all"\t"all"\t"$withPeak"\t"$denominator >> observedCounts.txt

species=rheMac10
denominator=$(cat ${species}_oCGIs_noFeatures.bed | wc -l)
species=rheMac2
for mark in ac me2
do
    bedtools intersect -wa -u -a ${species}_oCGIs_noFeatures.bed -b /home/ak2267/project/EnhancerClasses/${species}/${mark}/merge_*_overlap*_named.bed >> rheMac10_allIntersects.bed
done
mark=ac
bedtools intersect -wa -u -a ${species}_oCGIs_noFeatures.bed -b /home/ak2267/project/EnhancerClasses/Limb/${species}/merge_*_overlap_named.bed >> rheMac10_allIntersects.bed
species=rheMac10
withPeak=$(cut -f 4 ${species}_allIntersects.bed | sort -u | wc -l)
echo -e $species"\t"all"\t"all"\t"$withPeak"\t"$denominator >> observedCounts.txt

species=mm39
denominator=$(cat ${species}_oCGIs_noFeatures.bed | wc -l)
species=mm9
for mark in ac me2
do
    bedtools intersect -wa -u -a ${species}_oCGIs_noFeatures.bed -b /home/ak2267/project/EnhancerClasses/${species}/${mark}/merge_*_overlap*_named.bed >> mm39_allIntersects.bed
done
mark=ac
bedtools intersect -wa -u -a ${species}_oCGIs_noFeatures.bed -b /home/ak2267/project/EnhancerClasses/Limb/${species}/merge_*_overlap_named.bed >> mm39_allIntersects.bed
species=mm39
withPeak=$(cut -f 4 ${species}_allIntersects.bed | sort -u | wc -l)
echo -e $species"\t"all"\t"all"\t"$withPeak"\t"$denominator >> observedCounts.txt

rm *_allIntersects.bed


# run CGI-centric job file 2 (do shuffling)
cd /gpfs/gibbs/pi/noonan/ak2267/singleSpecies/shuffling

for species in rheMac10 calJac4 mm39 rn7 susScr11 canFam6 felCat9 equCab3 hg19; do mkdir intermediate_${species}; done
mkdir countsWithPeak

rm countsWithPeak/*
rm intermediate_*/*

# run rhesus as test
grep rheMac10 230112_step2_do_shuffling_jobFile.txt > 230112_step2_do_shuffling_jobFile_rheMac10.txt
dsq --job-file 230112_step2_do_shuffling_jobFile_rheMac10.txt --mem-per-cpu 10G -c 1 --mail-type FAIL,END --max-jobs 10
sbatch dsq-230112_step2_do_shuffling_jobFile_rheMac10-2023-01-13.sh # 19920841

# run rest of species
grep -v rheMac10 230112_step2_do_shuffling_jobFile.txt > 230116_step2_do_shuffling_jobFile_allRest.txt
dsq --job-file 230116_step2_do_shuffling_jobFile_allRest.txt --mem-per-cpu 5G -c 1 --mail-type FAIL,END --max-jobs 60
sbatch dsq-230116_step2_do_shuffling_jobFile_allRest-2023-01-16.sh # 19946961

# run all species at once on re-run
dsq --job-file 230112_step2_do_shuffling_jobFile.txt --mem-per-cpu 5G -c 1 --mail-type FAIL,END --max-jobs 100
sbatch dsq-230112_step2_do_shuffling_jobFile-2023-03-02.sh # 21576996
# ^^ this was a mistake because the rhesus ones fail with an out of memory error - will have to re-run rhesus

dsqa -j 21576996 -s OUT_OF_MEMORY -f 230112_step2_do_shuffling_jobFile.txt > 230305_step2_redo_jobFile.txt
dsq --job-file 230305_step2_redo_jobFile.txt --mem-per-cpu 10G -c 1 --mail-type FAIL,END --max-jobs 25
sbatch dsq-230305_step2_redo_jobFile-2023-03-05.sh # 21636920

# Six jobs timed out (I suspect they got stuck on slow nodes):
dsqa -j 21576996 -f 230112_step2_do_shuffling_jobFile.txt -s TIMEOUT
# susScr11_jobNum41
# susScr11_jobNum46
# susScr11_jobNum47
# canFam6_jobNum82
# canFam6_jobNum86
# canFam6_jobNum90

# move these output files to a special folder and re-run the jobs
mv countsWithPeak/susScr11_*_jobNum41.txt timedOutJobs_230312/
mv countsWithPeak/susScr11_*_jobNum46.txt timedOutJobs_230312/
mv countsWithPeak/susScr11_*_jobNum47.txt timedOutJobs_230312/
mv countsWithPeak/canFam6_*_jobNum82.txt timedOutJobs_230312/
mv countsWithPeak/canFam6_*_jobNum86.txt timedOutJobs_230312/
mv countsWithPeak/canFam6_*_jobNum90.txt timedOutJobs_230312/

rm intermediate_susScr11/*jobNum41*
rm intermediate_susScr11/*jobNum46*
rm intermediate_susScr11/*jobNum47*
rm intermediate_canFam6/*jobNum82*
rm intermediate_canFam6/*jobNum86*
rm intermediate_canFam6/*jobNum90*

dsqa -j 21576996 -f 230112_step2_do_shuffling_jobFile.txt -s TIMEOUT > 230312_step2_timedOutJobs_jobFile.txt
dsq --job-file 230312_step2_timedOutJobs_jobFile.txt --mem-per-cpu 10G -c 1 --mail-type FAIL,END
sbatch dsq-230312_step2_timedOutJobs_jobFile-2023-03-12.sh # 21707456

# summarize results using summarizeReshuffling_CGI-centric.py
for species in rheMac10 calJac4 mm39 rn7 susScr11 canFam6 felCat9 equCab3
do
    for tissue in brain liver muscle testis
    do
        for mark in H3K4me3 H3K27ac H3K4me1
        do
            python summarizeReshuffling_CGI-centric.py ${species} ${tissue} ${mark}
done; done; done >> reshufflingSummary_CGI-centric.txt

for species in hg19 rheMac2 mm9
do
    python summarizeReshuffling_CGI-centric.py ${species} devBrain ac
    python summarizeReshuffling_CGI-centric.py ${species} devBrain me2
    python summarizeReshuffling_CGI-centric.py ${species} devLimb ac
done >> reshufflingSummary_CGI-centric.txt

for species in hg19 rheMac10 calJac4 mm39 rn7 susScr11 canFam6 felCat9 equCab3
do
    python summarizeReshuffling_CGI-centric.py ${species} all all
done >> reshufflingSummary_CGI-centric.txt

# download reshufflingSummary_CGI-centric.txt for use in making R plots for Fig 1
cd /Users/acadiak/Desktop/CGI/singleSpecies
scp ak2267@ruddle.hpc.yale.edu:/gpfs/gibbs/pi/noonan/ak2267/singleSpecies/shuffling/reshufflingSummary_CGI-centric.txt .






