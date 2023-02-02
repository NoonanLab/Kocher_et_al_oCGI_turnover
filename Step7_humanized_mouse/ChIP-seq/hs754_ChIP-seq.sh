# 10/26/22
# ChIP-seq analysis from hs754 humanized mouse

# Purpose: map, call peaks, do differential peak calling, visualize ChIP data for hs754
# Includes samples from two separate ChIP-seq experiments (see hs754_ChIPseq_set1.sh and hs754_ChIPseq_set2.sh)
# Changes to this pipeline: remove multimapping and duplicate reads before peak calling

# Samples from set 1: H3K4me3 & H3K27me3 from e11.5 & e17.5, CTCF from e17.5
# Samples from set 2: H3K27ac from e11.5 & e17.5
# All samples are diencephalon

## SET 1
# data info is here: http://fcb.ycga.yale.edu:3010/ZtBQFftA3GLI9L8gi59dsmMAbVkK6/sample_dir_000007498/
# location on Ruddle: /ycga-gpfs/sequencers/illumina/sequencerC/runs/220215_A01519_0072_BH7HV3DSX3/Data/Intensities/BaseCalls/Unaligned-150/Project_Ak2267/
# see file Feb2022_ChIP_Kocher.xlsx for sample info

## SET 2
# data info is here: http://fcb.ycga.yale.edu:3010/roiZ960eWpD2w7SEB1qFtg3Ns0ph5/sample_dir_000007810/
# location on Ruddle: /ycga-gpfs/sequencers/illumina/sequencerD/runs/220328_A00124_0413_AHFMJFDSX3/Data/Intensities/BaseCalls/Unaligned-150/Project_Ak2267/
# see file Mar2022_ChIP_Kocher.xlsx for sample info

# work here
mkdir /home/ak2267/project/hs754_ChIP
cd /home/ak2267/project/hs754_ChIP

mkdir /home/ak2267/project/hs754_ChIP/fastqc
mkdir /home/ak2267/project/hs754_ChIP/counts

mkdir /home/ak2267/scratch60/hs754_ChIP
mkdir /home/ak2267/scratch60/hs754_ChIP/bam
mkdir /home/ak2267/scratch60/hs754_ChIP/bam/log
mkdir /home/ak2267/scratch60/hs754_ChIP/peaks
mkdir /home/ak2267/scratch60/hs754_ChIP/visualize
mkdir /home/ak2267/scratch60/hs754_ChIP/visualize/log

# MAKE NEW SAMPLE INFO FILE ENCOMPASSING EVERYTHING
cat /home/ak2267/project/hs754_ChIPseq_set1/hs754_ChIPseq_set1_sampleInfo.txt /home/ak2267/project/hs754_ChIPseq_set2/hs754_ChIPseq_set2_sampleInfo.txt > hs754_ChIP-seq_sampleInfo.txt

################### make humanized genome and bowtie indexes based on it (done previously in hs754_ChIPseq_set1.sh, here is code)
# make humanized genome - based in mm39
cd /home/ak2267/genomes
wget https://hgdownload.soe.ucsc.edu/goldenPath/mm39/bigZips/mm39.fa.gz

cd /home/ak2267/project/hs754_ChIPseq_set1
# code to get sequence, to put into makeHumanizedFasta.py

## /home/ak2267/project/hs754_ChIPseq_set1/bed/mouseCoord.bed
# mm39	chr13	72433164	72436073	5primeHomology	500	+
# mm39	chr13	72436073	72441214	humanizedMouse	500	+
# mm39	chr13	72441214	72444168	3primeHomology	500	+

## /home/ak2267/project/hs754_ChIPseq_set1/bed/humanCoord.bed
# hg38	chr5	3193946	3199477	humanizedHuman	500	-

bedtools getfasta -fi /home/ak2267/genomes/mm39.fa -bed bed/mouseCoord.bed
bedtools getfasta -s -fi /gpfs/ycga/datasets/genomes/Homo_sapiens/UCSC/hg38/Sequence/WholeGenomeFasta/genome.fa -bed bed/humanCoord.bed

# paste these sequences into makeHumanizedFasta.py

cd /home/ak2267/genomes
python /home/ak2267/project/hs754_ChIPseq_set1/makeHumanizedFasta.py mm39.fa mm39_humanizedHs754.fa

# check the replacement
cd /home/ak2267/project/hs754_ChIPseq_set1
bedtools getfasta -fi /home/ak2267/genomes/mm39_humanizedHs754.fa -bed bed/mouseCoord.bed

# make bowtie indexes
python /home/ak2267/project/hs754_ChIPseq_set1/prepFiles/prepFiles_bowtieIndexesHs754.py mm39
python /home/ak2267/project/hs754_ChIPseq_set1/prepFiles/prepFiles_bowtieIndexesHs754.py mm39_humanizedHs754
for i in $(ls *_bowtieIndexes.batch); do echo $i; sbatch $i; done
# 14166241, 14166242
####################


# CALCULATE EFFECTIVE GENOME SIZE for all genomes
# done already for mm39 in 220711_processRollerData_finalRun.sh
# https://deeptools.readthedocs.io/en/develop/content/feature/effectiveGenomeSize.html
# 'mm39':2309746861


# make job files with /home/ak2267/project/hs754_ChIP/prepFiles_hs754_ChIP-seq.py
cd /home/ak2267/project/hs754_ChIP
python prepFiles_hs754_ChIP-seq.py hs754_ChIP-seq_sampleInfo.txt

# makes the following:
# 221026_fastqc_jobFile.txt
# 221026_align_jobFile.txt
# 221026_sort_jobFile.txt
# 221026_filter_jobFile.txt
# 221026_callPeaks_jobFile.txt
# 221026_visualize_jobFile.txt
# 221026_count_jobFile.txt

# run FASTQC
dsq --job-file 221026_fastqc_jobFile.txt -J fastqc --mem-per-cpu 5G -c 1 --mail-type FAIL,END
sbatch dsq-221026_fastqc_jobFile-2022-10-26.sh # 18029565

rm fastqc/*.zip
tar -zcvf fastqc.gz fastqc/

# download to laptop to view html files
cd /Users/acadiak/Desktop/hs754-ChIP/Sequencing/fastqc
scp ak2267@ruddle.hpc.yale.edu:/home/ak2267/project/hs754_ChIP/fastqc.gz .

# run ALIGN
dsq --job-file 221026_align_jobFile.txt -J align --mem-per-cpu 5G -c 8 --mail-type FAIL,END
sbatch dsq-221026_align_jobFile-2022-10-26.sh # 18029883

# run SORT
dsq --job-file 221026_sort_jobFile.txt -J sort --mem-per-cpu 5G -c 2 --mail-type FAIL,END
sbatch dsq-221026_sort_jobFile-2022-10-28.sh # 18098246

# run FILTER
dsq --job-file 221026_filter_jobFile.txt -J filter --mem-per-cpu 5G -c 1 --mail-type FAIL,END
sbatch dsq-221026_filter_jobFile-2022-10-28.sh # 18098489

# remove unneeded intermediate bam files (keep pre-sort and pre-filter for now)
rm ~/scratch60/hs754_ChIP/bam/*_intermediateFiltered.bam
rm ~/scratch60/hs754_ChIP/bam/*_sorted.bam

# CALL PEAKS
dsq --job-file 221026_callPeaks_jobFile.txt -J callPeaks --mem-per-cpu 5G -c 1 --mail-type FAIL,END
sbatch dsq-221026_callPeaks_jobFile-2022-10-28.sh # 18098573

for i in $(ls *narrowPeak)
do
    echo ${i}
    bedtools intersect -wa -a ${i} -b /home/ak2267/project/hs754_ChIPseq_set1/bed/mouseCoord.bed
done

# VISUALIZE (make bigWigs and bigBeds)
dsq --job-file 221026_visualize_jobFile.txt -J visualize --mem-per-cpu 5G -c 2 --mail-type FAIL,END
sbatch dsq-221026_visualize_jobFile-2022-10-28.sh # 18098615

# move bigwigs to server for viewing on the browser
cd /home/ak2267/scratch60/hs754_ChIP/visualize
for i in $(ls *.bw)
do
    scp ${i} ak2267@10.5.37.220:/home/ak2267/akocher_www/hs754_ChIP/mm39/
    sleep 5s
done
ssh ak2267@10.5.37.220

# move bigBeds to the browser
cd /home/ak2267/scratch60/hs754_ChIP/visualize
for i in $(ls *.bb)
do
    scp ${i} ak2267@10.5.37.220:/home/ak2267/akocher_www/hs754_ChIP/mm39/
    sleep 1s
done

# also transfer mouseCoord_mm39.bed
# mm39	chr13	72433164	72436073	5primeHomology	500	+
# mm39	chr13	72436073	72441214	humanizedRegion	500	+
# mm39	chr13	72441214	72444168	3primeHomology	500	+
# and humanCoord_mm39.bed
# mm39	chr13	72436073	72441214	humanizedCoordinates	500	+
bedToBigBed mouseCoord_mm39.bed /home/ak2267/genomes/chrom.sizes/mm39.chrom.sizes mouseCoord_mm39.bb
scp mouseCoord_mm39.bb ak2267@10.5.37.220:/home/ak2267/akocher_www/hs754_ChIP/mm39/

bedToBigBed humanCoord_mm39.bed /home/ak2267/genomes/chrom.sizes/mm39.chrom.sizes humanCoord_mm39.bb
scp humanCoord_mm39.bb ak2267@10.5.37.220:/home/ak2267/akocher_www/hs754_ChIP/mm39/

# make trackDb.txt file for track hub
cd /home/ak2267/project/hs754_ChIP
python makeTrackHubFiles.py > trackDb.txt
scp trackDb.txt ak2267@10.5.37.220:/home/ak2267/akocher_www/hs754_ChIP/mm39/

### GENERATE UNION OF PEAKS BETWEEN GENOTYPES

# The issue: everything downstream of hs754 is offset by 390bp, and peaks within hs754 are also not ideally aligned (e.g. CTCF)
# Adjusting peaks that are fully downstream is easy - just need to be shifted back 390bp. Strategy will be to shift them back, merge in mm39, and then shift forward for HTSeq in HUM
# More complicated:
# 1) peaks that start before hs754 and end within it
# 2) peaks that start before hs754 and end downstream of it by > 390bp
# 2.5) peaks that start before hs54 and end downstream but by less than 390bp)
# 3) peaks contained fully within hs754 (i.e. 5531bp after the end of the mouse 5' homology arm)
# 4) peaks that start within hs754 and end downstream of it
# 5) peaks that start within 390bp downstream of hs754 - should they be shifted back the full 390bp?
# 6) total peaks on chr13

# transfer contents of /Users/acadiak/Desktop/hs754-ChIP/Sequencing/PeakShifts to Ruddle
cd /Users/acadiak/Desktop/hs754-ChIP/Sequencing/PeakShifts
scp summarizePeakCategories.py ak2267@ruddle.hpc.yale.edu:/home/ak2267/project/hs754_ChIP/peakShifts
scp moveMm39_toHg38.py ak2267@ruddle.hpc.yale.edu:/home/ak2267/project/hs754_ChIP/peakShifts
scp adjustHUMpeaksInMm39.py ak2267@ruddle.hpc.yale.edu:/home/ak2267/project/hs754_ChIP/peakShifts
scp makeGTF_newPeakNumbers.py ak2267@ruddle.hpc.yale.edu:/home/ak2267/project/hs754_ChIP/
scp moveHg38_toHUMmm39.py ak2267@ruddle.hpc.yale.edu:/home/ak2267/project/hs754_ChIP/peakShifts
scp adjustMergedPeaks_inHUMmm39.py ak2267@ruddle.hpc.yale.edu:/home/ak2267/project/hs754_ChIP/peakShifts

# First step: identify the number of peaks in each of these categories
rm /home/ak2267/scratch60/hs754_ChIP/peaks/*.gappedPeak
for i in $(ls /home/ak2267/scratch60/hs754_ChIP/peaks/*HUM*_peaks.*Peak)
do
    python peakShifts/summarizePeakCategories.py ${i}
done

# Conclusion: there are only humanized peaks in categories 2 (H3K27me3), 3 (H3K27ac, H3K4me3, CTCF), 4 (H3K27ac at e17.5), and 5 (H3K27ac at e17.5)

# HOWEVER if I merge peaks within WT & HUM and then re-assess...

for genotype in WT HUM
do
    # MULTIPLEX A
    multiplex=A
    for mark in H3K4me3 H3K27me3 CTCF
    do
        cat /home/ak2267/scratch60/hs754_ChIP/peaks/${multiplex}_e17.5_${mark}_${genotype}_*_peaks.*Peak > mergedPeaks/${multiplex}_${mark}_${genotype}_cat.bed
        sort -k1,1 -k2,2n mergedPeaks/${multiplex}_${mark}_${genotype}_cat.bed > mergedPeaks/${multiplex}_${mark}_${genotype}_sort.bed
        bedtools merge -i mergedPeaks/${multiplex}_${mark}_${genotype}_sort.bed -c 4 -o collapse > mergedPeaks/${multiplex}_e17.5_${mark}_${genotype}_merged.bed
    done
    # MULTIPLEX B
    multiplex=B; mark=H3K4me3
    cat /home/ak2267/scratch60/hs754_ChIP/peaks/${multiplex}_e11.5_${mark}_${genotype}_*_peaks.*Peak > mergedPeaks/${multiplex}_${mark}_${genotype}_cat.bed
    sort -k1,1 -k2,2n mergedPeaks/${multiplex}_${mark}_${genotype}_cat.bed > mergedPeaks/${multiplex}_${mark}_${genotype}_sort.bed
    bedtools merge -i mergedPeaks/${multiplex}_${mark}_${genotype}_sort.bed -c 4 -o collapse > mergedPeaks/${multiplex}_e11.5_${mark}_${genotype}_merged.bed
    # MULTIPLEX C
    multiplex=C; mark=H3K27me3
    cat /home/ak2267/scratch60/hs754_ChIP/peaks/${multiplex}_e11.5_${mark}_${genotype}_*_peaks.*Peak > mergedPeaks/${multiplex}_${mark}_${genotype}_cat.bed
    sort -k1,1 -k2,2n mergedPeaks/${multiplex}_${mark}_${genotype}_cat.bed > mergedPeaks/${multiplex}_${mark}_${genotype}_sort.bed
    bedtools merge -i mergedPeaks/${multiplex}_${mark}_${genotype}_sort.bed -c 4 -o collapse > mergedPeaks/${multiplex}_e11.5_${mark}_${genotype}_merged.bed
    # MULTIPLEX D
    multiplex=D; mark=H3K27ac
    cat /home/ak2267/scratch60/hs754_ChIP/peaks/${multiplex}_e17.5_${mark}_${genotype}_*_peaks.*Peak > mergedPeaks/${multiplex}_${mark}_${genotype}_cat.bed
    sort -k1,1 -k2,2n mergedPeaks/${multiplex}_${mark}_${genotype}_cat.bed > mergedPeaks/${multiplex}_${mark}_${genotype}_sort.bed
    bedtools merge -i mergedPeaks/${multiplex}_${mark}_${genotype}_sort.bed -c 4 -o collapse > mergedPeaks/${multiplex}_e17.5_${mark}_${genotype}_merged.bed
    # MULTIPLEX E
    multiplex=E; mark=H3K27ac
    cat /home/ak2267/scratch60/hs754_ChIP/peaks/${multiplex}_e11.5_${mark}_${genotype}_*_peaks.*Peak > mergedPeaks/${multiplex}_${mark}_${genotype}_cat.bed
    sort -k1,1 -k2,2n mergedPeaks/${multiplex}_${mark}_${genotype}_cat.bed > mergedPeaks/${multiplex}_${mark}_${genotype}_sort.bed
    bedtools merge -i mergedPeaks/${multiplex}_${mark}_${genotype}_sort.bed -c 4 -o collapse > mergedPeaks/${multiplex}_e11.5_${mark}_${genotype}_merged.bed
done

for i in $(ls mergedPeaks/*HUM*_merged.bed)
do
    python peakShifts/summarizePeakCategories.py ${i}
done

# eliminates site in category 5 - so now just categories 2 (H3K27me3), 3 (H3K27ac, H3K4me3, CTCF), 4 (H3K27ac at e17.5)

# Strategy for category 2: shift end of peak back 390bp and that's the coordinates in mm39
# Strategy for category 3: assign start & end coordinates in hg38 and use liftOver to lift to mm39
# Strategy for category 4: assign start coordinates in hg38 and use liftOver to lift to mm39; subtract 390 from peak end

# Steps in this updated pipeline for differential peak assessment:
# 1) call peaks in each sample, in either WT or HUM coordinates
# 2) merge within each genotype
# 3) shift HUM peaks back to WT coordinate system in order to make a merged set across genotypes

#    A) peaks that start upstream and end downstream: subtract 390bp from end
#    B) peaks contained within humanized region: put into hg38 coordinates, liftOver to mm39
#    C) peaks that start within and end downstream: liftOver start to mm39, subtract 390bp from end
#    D) peaks downstream of humanized region: subtract 390bp from both start and end

# 4) THEN count reads in merged peaks (which are different coordinates in the WT vs HUM) and proceed with DESeq2


# first do the liftOver step (extract peaks in hum region & put in hg38 coordinates, lift to mm39 and put in bed file)
# only extracts peaks that are contained fully within the humanized coordinates (i.e. nothing for K27me3 peaks)
for i in $(ls mergedPeaks/*HUM*_merged.bed)
do
    # extract peaks falling completely within hs754 using bedtools intersect (i.e. cat 3 but not 2 or 4)
    sampleName=$(basename ${i}} | awk -F'[_]' '{print $1"_"$2"_"$3"_"$4}')
    echo ${sampleName}
    bedtools intersect -f 1 -wa -a ${i} -b humanCoord_mm39.bed > peakShifts/${sampleName}_humanizedInMm39.bed
    
    # put these peaks into hg38 coordinates by reversing them
    # also remove last few columns of bed so it's only chr start end name
    python peakShifts/moveMm39_toHg38.py peakShifts/${sampleName}_humanizedInMm39.bed > peakShifts/${sampleName}_inHg38.bed
    
    # liftOver to mm39
    liftOver -minMatch=0.3 peakShifts/${sampleName}_inHg38.bed /home/ak2267/genomes/chain/hg38ToMm39.over.chain.gz peakShifts/${sampleName}_hg38LOmm39.bed unMapped
done

# CATEGORY 4 in file mergedPeaks/D_e17.5_H3K27ac_HUM_merged.bed - do this manually
# chr13	72440477	72441604    D_e17.5_H3K27ac_HUM_3_peak_12406,D_e17.5_H3K27ac_HUM_2_peak_17050,D_e17.5_H3K27ac_HUM_1_peak_16832,D_e17.5_H3K27ac_HUM_2_peak_17051

# coord for peak that hangs over:           chr13	72440477	72441798
# coord ending with end of HUM region:      chr13	72440477    72441604
# coord in hg38 (flip and go from end):     chr5    3193946 3195073
#   startInHumanCoord = humanRegionStart + (mouseRegionEnd + 390 - humanizedMousePeakEnd) = humanRegionStart (3193946)
#   endInHumanCoord = humanRegionEnd (3199477) - (humanizedMousePeakStart (72440477) - mouseRegionStart (72436073))
# coord after lifting from hg38 to mm39:    chr13	72440000	72441044
# FINAL COORD                               chr13	72440000	72441408
# after using this 5' end but using original 3' end minus 390bp (72441798 - 390)

# manually edit e17.5_H3K27ac_HUM_2_hg38LOmm39.bed to add this entry
# chr13	72440000	72441408    D_e17.5_H3K27ac_HUM_3_peak_12406,D_e17.5_H3K27ac_HUM_2_peak_17050,D_e17.5_H3K27ac_HUM_1_peak_16832,D_e17.5_H3K27ac_HUM_2_peak_17051

# python script to integrate them and also renumber peaks in other categories
for i in $(ls mergedPeaks/*HUM*_merged.bed)
do
    sampleName=$(basename ${i}} | awk -F'[_]' '{print $1"_"$2"_"$3"_"$4}')
    echo ${sampleName}
    python peakShifts/adjustHUMpeaksInMm39.py ${i} peakShifts/${sampleName}_hg38LOmm39.bed > peakShifts/${sampleName}_HUMpeaksInMm39.bed
done

# 3) make merged set in WT coordinates with peaks named for coordinates in WT (will make things easier in DESeq2)
for i in A_e17.5_CTCF A_e17.5_H3K27me3 A_e17.5_H3K4me3 B_e11.5_H3K4me3 C_e11.5_H3K27me3 D_e17.5_H3K27ac E_e11.5_H3K27ac
do
    echo $i
    cat mergedPeaks/${i}_WT_merged.bed peakShifts/${i}_HUM_HUMpeaksInMm39.bed > mergedInMm39/${i}_cat.bed
    cut -f 1,2,3,4 mergedInMm39/${i}_cat.bed | sort -k1,1 -k2,2n - > mergedInMm39/${i}_sort.bed
    bedtools merge -o collapse -c 4 -i mergedInMm39/${i}_sort.bed > mergedInMm39/${i}_mergedPeaks.bed
    python makeGTF_newPeakNumbers.py mergedInMm39/${i}_mergedPeaks.bed mergedInMm39/${i}_mergedPeaks.gtf mergedInMm39/${i}_mergedPeaks_renamed.bed
done

# 4) shift the merged set (in GTF form) back to HUM coordinates by doing the inverse of what you did above in step 2)
for i in $(ls mergedInMm39/*_mergedPeaks_renamed.bed)
do
    # extract peaks falling completely within hs754 using bedtools intersect
    sampleName=$(basename ${i}} | awk -F'[_]' '{print $1"_"$2"_"$3}')
    echo ${sampleName}
    bedtools intersect -f 1 -wa -a ${i} -b humanCoord_mm39.bed > peakShifts/${sampleName}_mergedPeaks_renamed_fallingInHs754.bed

    # liftOver to mm39
    liftOver -minMatch=0.3 peakShifts/${sampleName}_mergedPeaks_renamed_fallingInHs754.bed /home/ak2267/genomes/chain/mm39ToHg38.over.chain.gz peakShifts/${sampleName}_mergedPeaks_renamed_fallingInHs754_LOhg38.bed unMapped
    
    # put these peaks into humanized mm39 coordinates by reversing them
    python peakShifts/moveHg38_toHUMmm39.py peakShifts/${sampleName}_mergedPeaks_renamed_fallingInHs754_LOhg38.bed > peakShifts/${sampleName}_mergedPeaks_renamed_fallingInHs754_LOhg38_backToHUM.bed
done

# CATEGORY 4 in file mergedPeaks/D_e17.5_H3K27ac_HUM_merged.bed - do this manually
bedtools intersect -a mergedInMm39/D*_mergedPeaks_renamed.bed -b humanCoord_mm39.bed
# coord in WT mouse:                               chr13	72440000	72441753	chr13:72440000-72441753
# coord ending with end of humanized region:       chr13	72440000	72441214
# coord after lifting to hg38:                     chr5	3193901	3195052
# coord in HUMmm39:                                chr13	72440498	72441994
#    startInHUMmouseCoord = mouseRegionStart (72436073) + (humanRegionEnd (3199477) - peakInHumanEnd (3195052))
#    endInHUMmouseCoord = simply end from WT coord + 390
# FINAL COORD:                                     chr13	72440477	72441994    chr13:72440000-72441753
# adjust start to be 72440477 like 5'-most peak in HUM coordinates (liftOver shifted it by 20bp)

# manually add this entry to peakShifts/D_e17.5_H3K27ac_mergedPeaks_renamed_fallingInHs754_LOhg38_backToHUM.bed
# 

# python script to integrate them and also renumber peaks in other categories, with 'P###" names in col 4 to easily make a GTF file for HTSeq
for i in A_e17.5_CTCF A_e17.5_H3K27me3 A_e17.5_H3K4me3 B_e11.5_H3K4me3 C_e11.5_H3K27me3 D_e17.5_H3K27ac E_e11.5_H3K27ac
do
    echo ${i}
    python peakShifts/adjustMergedPeaks_inHUMmm39.py mergedInMm39/${i}_mergedPeaks_renamed.bed peakShifts/${i}_mergedPeaks_renamed_fallingInHs754_LOhg38_backToHUM.bed > mergedInHUMmm39/${i}_mergedPeaks_inHUMmm39.bed
    python ~/Scripts/makeGTF_nameFromCol4.py mergedInHUMmm39/${i}_mergedPeaks_inHUMmm39.bed mergedInHUMmm39/${i}_mergedPeaks_inHUMmm39.gtf
done

# FILES FOR HTSeq:
# WT:   mergedInMm39/${multiplexGroup}_${timePoint}_${mark}_mergedPeaks.gtf
# HUM:  mergedInHUMmm39/${multiplexGroup}_${timePoint}_${mark}_mergedPeaks_inHUMmm39.gtf 

# run COUNT with HTSeq
dsq --job-file 221026_count_jobFile.txt -J count --mem-per-cpu 5G -c 1 --mail-type FAIL,END
sbatch dsq-221026_count_jobFile-2022-10-30.sh # 18126404

# run COUNT with HTSeq - this time with -r pos (reads are sorted by position, not name - name is the default so it thought all my paired reads lacked a mate)
dsq --job-file 221026_count_jobFile.txt -J count --mem-per-cpu 5G -c 1 --mail-type FAIL,END
sbatch dsq-221026_count_jobFile-2022-10-31.sh # 18129073
# THIS FAILED BECAUSE OF LOW MEMORY

# NAME SORT filtered bam files 
dsq --job-file 221026_nameSort_jobFile.txt -J count --mem-per-cpu 5G -c 1 --mail-type FAIL,END
sbatch dsq-221026_nameSort_jobFile-2022-11-02.sh # 18168810

# run COUNT - for each mark, also run on input sequence
# results are in /home/ak2267/project/hs754_ChIP/counts with format:
# A_e17.5_H3K4me3_WT_1.quant, and its input, A_e17.5_inputH3K4me3_WT_1.quant 
dsq --job-file 221026_count_jobFile.txt -J count --mem-per-cpu 5G -c 1 --mail-type FAIL,END
sbatch dsq-221026_count_jobFile-2022-11-02.sh # 18169163

# this worked much better! more counts per peak. I think before it was throwing out a lot of reads because it couldn't find the pair

# gzip and download HTSeq files for use with DESeq2 in R (hs754_ChIP-seq.R)
cd /home/ak2267/project/hs754_ChIP
tar -zcvf counts_hs754_ChIP.gz counts
/Users/acadiak/Desktop/CGI/hs754/ChIP
scp ak2267@ruddle.hpc.yale.edu:/home/ak2267/project/hs754_ChIP/counts_hs754_ChIP.gz .
# HTSeq files are in /home/ak2267/project/hs754_ChIP/counts




