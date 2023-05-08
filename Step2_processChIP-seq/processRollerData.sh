# 7/11/22
# Purpose: process all Roller, Stamper et al ChIP data
# download, align, call peaks, identify replicating peaks

# work here
cd /home/ak2267/project/Roller/ChIP/batch

##################
# submit batch files that DOWNLOAD genome fastas (unmasked) and MAKE BOWTIE2 INDEXES in /gpfs/gibbs/pi/noonan/ak2267/genomes
# this takes 20-30 min each

for species in rheMac10 calJac4 mm39 rn7 oryCun2 susScr11 canFam6 felCat9 equCab3
do
    python prepFiles/prepFiles_bowtieIndexes_finalRun.py ${species}
done
cd /home/ak2267/project/Roller/ChIP/batch

mv *_bowtieIndexes.batch bowtieIndexes/
for i in $(ls bowtieIndexes/*_bowtieIndexes.batch); do echo $i; sbatch $i; done
# Submitted batch job 15610892,15610893,15610894,15610895,15610896,15610897,15610898,15610899,15610900

# CALCULATE EFFECTIVE GENOME SIZE for all genomes
# using unique-kmers.py
# https://deeptools.readthedocs.io/en/develop/content/feature/effectiveGenomeSize.html

# install khmer
ml miniconda
mamba create -yn khmerEnv python khmer
conda activate khmerEnv
# then can run program unique-kmers.py

ml miniconda
conda activate khmerEnv
# test on rheMac10
unique-kmers.py -k 50 /gpfs/gibbs/pi/noonan/ak2267/genomes/

# run unique-kmers.py
for species in rheMac10 calJac4 mm39 rn7 susScr11 canFam6 felCat9 equCab3
do
    echo $species
    unique-kmers.py -k 50 /gpfs/gibbs/pi/noonan/ak2267/genomes/${species}.fa.gz
done

# store results in this dictionary - use in python script above that generates bowtieIndex job file (prepFiles_bowtieIndexes_finalRun.py)
# Effective_GenomeSize = {'calJac4':2621714877,'canFam6':2240314204,'equCab3':2365156725,'felCat9':2352800598,'rheMac10':2653677440,'mm39':2309746861,'rn7':2369256514,'susScr11':2359453421}

# gunzip the genomes for use with getfasta / faCount later on in the pipeline
for species in rheMac10 calJac4 mm39 rn7 susScr11 canFam6 felCat9 equCab3
do
    echo $species
    gunzip ${species}.fa.gz
done

##################
# MAKE JOB FILES FOR DOWNLOAD, ALIGN, CALL PEAKS
cd /home/ak2267/project/Roller/ChIP/batch
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-7127/E-MTAB-7127.sdrf.txt

module load dSQ
python prepFiles/prepFiles_processRoller_finalRun.py E-MTAB-7127.sdrf.txt
# makes
220712_downloadFastq_jobFile.txt
220712_fastqc_jobFile.txt
220712_align_jobFile.txt
220713_sort_jobFile.txt
220713_filter_jobFile.txt
220713_callPeaks_jobFile.txt

# DOWNLOAD to /gpfs/gibbs/pi/noonan/ak2267/Roller/fastq
dsq --job-file 220712_downloadFastq_jobFile.txt --mem-per-cpu 5G -c 1 --mail-type FAIL,END
sbatch dsq-220712_downloadFastq_jobFile-2022-07-12.sh # 15612960
sbatch dsq-220712_downloadFastq_jobFile-2023-02-11.sh # 21094756 - Feb 2023

# run fastqc on fastq files
dsq --job-file 220712_fastqc_jobFile.txt --mem-per-cpu 5G -c 1 --mail-type FAIL,END
sbatch dsq-220712_fastqc_jobFile-2022-07-12.sh # 15614407

cd /home/ak2267/project/Roller/ChIP/batch
# ALIGN to /gpfs/gibbs/pi/noonan/ak2267/Roller/bam - keep unmapped reads
grep 'H3K4me3' 220712_align_jobFile.txt > 220712_alignH3K4me3_jobFile.txt
grep 'input' 220712_align_jobFile.txt > 220712_alignInput_jobFile.txt
grep 'H3K27ac' 220712_align_jobFile.txt > 220712_alignH3K27ac_jobFile.txt
grep 'H3K4me1' 220712_align_jobFile.txt > 220712_alignH3K4me1_jobFile.txt

# submit jobs
dsq --job-file 220712_alignH3K4me3_jobFile.txt --mem-per-cpu 5G -c 8 --mail-type FAIL,END
sbatch dsq-220712_alignH3K4me3_jobFile-2022-07-12.sh # 15615056
sbatch dsq-220712_alignH3K4me3_jobFile-2023-02-11.sh # 21095422 - Feb 2023
dsq --job-file 220712_alignInput_jobFile.txt --mem-per-cpu 5G -c 8 --mail-type FAIL,END
sbatch dsq-220712_alignInput_jobFile-2022-07-12.sh # 15615260
sbatch dsq-220712_alignInput_jobFile-2023-02-11.sh # 21096231 - Feb 2023
dsq --job-file 220712_alignH3K27ac_jobFile.txt --mem-per-cpu 5G -c 8 --mail-type FAIL,END
sbatch dsq-220712_alignH3K27ac_jobFile-2022-07-12.sh # 15615262
sbatch dsq-220712_alignH3K27ac_jobFile-2023-02-12.sh # 21102057
dsq --job-file 220712_alignH3K4me1_jobFile.txt --mem-per-cpu 5G -c 8 --mail-type FAIL,END
sbatch dsq-220712_alignH3K4me1_jobFile-2022-07-12.sh # 15615264
sbatch dsq-220712_alignH3K4me1_jobFile-2023-02-12.sh # 21099800

# DELETE FASTQ
rm /gpfs/gibbs/pi/noonan/ak2267/Roller/fastq/*H3K4me3*.fastq.gz
rm /gpfs/gibbs/pi/noonan/ak2267/Roller/fastq/*input*.fastq.gz
rm /gpfs/gibbs/pi/noonan/ak2267/Roller/fastq/*H3K27ac*.fastq.gz
rm /gpfs/gibbs/pi/noonan/ak2267/Roller/fastq/*H3K4me1*.fastq.gz

# SORT BAM FILES & DELETE UNSORTED
cd /home/ak2267/project/Roller/ChIP/batch
dsq --job-file 220713_sort_jobFile.txt --mem-per-cpu 5G -c 2 --mail-type FAIL,END
sbatch dsq-220713_sort_jobFile-2022-07-13.sh # 15616589
sbatch dsq-220713_sort_jobFile-2023-02-13.sh # 21105259
rm /gpfs/gibbs/pi/noonan/ak2267/Roller/bam/*_bowtie2.bam

# REMOVE MULTI-MAPPING READS AND DUPLICATES USING SAMBAMBA
# https://hbctraining.github.io/Intro-to-ChIPseq/lessons/03_align_and_filtering.html
dsq --job-file 220713_filter_jobFile.txt --mem-per-cpu 5G -c 2 --mail-type FAIL,END
sbatch dsq-220713_filter_jobFile-2022-07-13.sh # 15617124
sbatch dsq-220713_filter_jobFile-2023-02-13.sh # 21108268

# remove intermediate bam files
rm /gpfs/gibbs/pi/noonan/ak2267/Roller/bam/*_intermediateFiltered.bam

# COMBINE SEPARATED INPUT FILES
cd /home/ak2267/project/Roller/ChIP/batch/combineInputs
python ../prepFiles/prepFiles_combineInputs_finalRun.py
for i in $(ls *.batch); do echo $i; sbatch $i; done
# Submitted batch job 15617773,15617774,15617775,15617776,15617777,15617778,15617779,15617780,15617781,15617782
# 6 completed, 4 failed (all oryCun that I forgot to take out)
# Submitted batch job 21109470,21109471,21109472,21109473,21109474,21109475

# INDEX NEW INPUT FILES
cd /gpfs/gibbs/pi/noonan/ak2267/Roller/bam/
for file in canFam6_muscle_input_3 equCab3_testis_input_2 felCat9_muscle_input_2 rheMac10_muscle_input_2 rheMac10_muscle_input_3 susScr11_muscle_input_1
do
    echo $file
    samtools index ${file}_bowtie2_filtered.bam
done

# run Multiqc to get reports
cd /gpfs/gibbs/pi/noonan/ak2267/Roller/bam/
module load MultiQC/1.10.1-foss-2020b-Python-3.8.6
multiqc .
# download for viewing
cd /Users/acadiak/Desktop/H3K4me3/Roller/Apr2022/MultiQC

# CALL PEAKS
# some may fail due to too many chromosomes and will need a special flag in the command
# --buffer-size 1000 (default is 100000)
# see discussion here https://github.com/macs3-project/MACS/issues/313
cd /home/ak2267/project/Roller/ChIP/batch/
grep 'H3K4me3' 220713_callPeaks_jobFile.txt > 220713_callPeaksH3K4me3_jobFile.txt
grep 'H3K27ac' 220713_callPeaks_jobFile.txt > 220713_callPeaksH3K27ac_jobFile.txt
grep 'H3K4me1' 220713_callPeaks_jobFile.txt > 220713_callPeaksH3K4me1_jobFile.txt

dsq --job-file 220713_callPeaksH3K4me3_jobFile.txt --mem-per-cpu 10G -c 1 --mail-type FAIL,END
sbatch dsq-220713_callPeaksH3K4me3_jobFile-2022-07-13.sh # 15618399

dsq --job-file 220713_callPeaksH3K27ac_jobFile.txt --mem-per-cpu 10G -c 1 --mail-type FAIL,END
sbatch dsq-220713_callPeaksH3K27ac_jobFile-2022-07-13.sh # 15618277

dsq --job-file 220713_callPeaksH3K4me1_jobFile.txt --mem-per-cpu 10G -c 1 --mail-type FAIL,END
sbatch dsq-220713_callPeaksH3K4me1_jobFile-2022-07-13.sh # 15618154

# VISUALIZE - files going into /gpfs/gibbs/pi/noonan/ak2267/Roller/visualize
cd /home/ak2267/project/Roller/ChIP/batch/
grep 'H3K4me3' 220714_visualize_jobFile.txt > 220714_visualizeH3K4me3_jobFile.txt
grep 'H3K27ac' 220714_visualize_jobFile.txt > 220714_visualizeH3K27ac_jobFile.txt
grep 'H3K4me1' 220714_visualize_jobFile.txt > 220714_visualizeH3K4me1_jobFile.txt
grep 'input' 220714_visualize_jobFile.txt > 220714_visualizeInput_jobFile.txt

dsq --job-file 220714_visualizeH3K4me3_jobFile.txt --mem-per-cpu 5G -c 2 --mail-type FAIL,END
sbatch dsq-220714_visualizeH3K4me3_jobFile-2022-07-14.sh # 15626356
dsq --job-file 220714_visualizeH3K27ac_jobFile.txt --mem-per-cpu 5G -c 2 --mail-type FAIL,END
sbatch dsq-220714_visualizeH3K27ac_jobFile-2022-07-14.sh # 15626664
dsq --job-file 220714_visualizeH3K4me1_jobFile.txt --mem-per-cpu 5G -c 2 --mail-type FAIL,END
sbatch dsq-220714_visualizeH3K4me1_jobFile-2022-07-14.sh # 15626740
dsq --job-file 220714_visualizeInput_jobFile.txt --mem-per-cpu 5G -c 2 --mail-type FAIL,END
sbatch dsq-220714_visualizeInput_jobFile-2022-07-24.sh # 15676096

rm /gpfs/gibbs/pi/noonan/ak2267/Roller/visualize/temp_*.bed

##### Move files to the server

# Transfer genomes.txt and trackDb_{species}.txt files from ruddle to the server
# df -h .
cd /home/ak2267/project/Roller/ChIP/batch
scp genomes.txt ak2267@10.5.37.220:/home/ak2267/akocher_www/CGI
for species in rheMac10 calJac4 mm39 rn7 susScr11 canFam6 felCat9 equCab3
do
    scp trackDb_${species}.txt ak2267@10.5.37.220:/home/ak2267/akocher_www/CGI/${species}
done

cd /gpfs/gibbs/pi/noonan/ak2267/Roller/visualize/

# have to cycle through each species & tissue individually because if I try to run more files I get an error message after ~30 files
#for species in rheMac10 calJac4 mm39 rn7 susScr11 canFam6 felCat9 equCab3
#for tissue in brain liver muscle testis

# swap out each mark in turn since jobs were still running when I started transferring - H3K4me3 H3K27ac H3K4me1 input
for species in rheMac10 calJac4 mm39 rn7 susScr11 canFam6 felCat9 equCab3
do
    for tissue in brain liver muscle testis
    do
	for i in $(ls ${species}*${tissue}*input*)
	do
	    scp ${i} ak2267@10.5.37.220:/home/ak2267/akocher_www/CGI/${species}
	done
	sleep 30s
    done;done


###############
# 7/27/22
# Count peaks in each file to see outliers - analysis in DataQuality_Roller.xlsx
cd /gpfs/gibbs/pi/noonan/ak2267/Roller/peaks/

for species in rheMac10 calJac4 mm39 rn7 susScr11 canFam6 felCat9 equCab3
do
    for tissue in brain liver muscle testis
    do
	for mark in H3K4me3 H3K27ac H3K4me1
	do
	    rep1=$(cat ${species}_${tissue}_${mark}_1_peaks.*Peak | wc -l)
	    rep2=$(cat ${species}_${tissue}_${mark}_2_peaks.*Peak | wc -l)
	    rep3=$(cat ${species}_${tissue}_${mark}_3_peaks.*Peak | wc -l)
	    echo $species $tissue $mark $rep1 $rep2 $rep3
	done >> peakSummary.txt
    done
done

#### Generate intersection of all replicates, doing a quality control step to remove outliers
# 9/13/22
# Remove a replicate if:
# it has >=50% MORE peaks than average of other two replicates
# OR
# if it has >= 20% FEWER peaks than average of other two replicates

# also rename col 4 with the peak name: species_tissue_mark_line#
# script outputs repsToUseRoller.txt (simply keeps track of which reps are good, for use downstream) and intersectCommands.txt (for batch submission)

python makeIntersectCommands.py peakSummary.txt repsToUseRoller.txt > intersectCommands.txt
dsq --job-file intersectCommands.txt --mem-per-cpu 5G -c 1 --mail-type FAIL,END
sbatch dsq-intersectCommands-2022-09-20.sh # 17135377
# results are here:
# /gpfs/gibbs/pi/noonan/ak2267/Roller/peaks/intersection/${species}_${tissue}_${mark}_intersection.bed


