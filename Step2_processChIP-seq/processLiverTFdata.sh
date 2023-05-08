# 7/20/22
# Purpose: download & process liver TF data 

# Data from Schmidt et al 2010 Science: CEBPA and HNF4a in liver from 5 animals: Human, Mouse, Dog, Opossum, Chicken
# Note human sample is hepatocytes not whole liver

# FOCUS ON THIS, BALLESTER 2014
# Data from Ballester et al 2014 eLife: CEBPA, HNF4a, ONECUT1, FOXA1 in liver from 5 mammals: Human, Rhesus, Dog, Mouse, Rat - some human samples are from older papers and therefore are hepatocytes (CEBPA, HNF4a)
# Focus on Rhesus, Dog, Mouse, Rat to match Roller 2021 genomes
# https://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-1509/

# work here
mkdir /home/ak2267/project/LiverTF/
mkdir /home/ak2267/project/LiverTF/batch
mkdir /home/ak2267/project/LiverTF/batch/prepFiles
mkdir /home/ak2267/project/LiverTF/batch/dsq
cd /home/ak2267/project/LiverTF/batch

# and put files into here
mkdir /gpfs/gibbs/pi/noonan/ak2267/LiverTF
mkdir /gpfs/gibbs/pi/noonan/ak2267/LiverTF/fastq
mkdir /gpfs/gibbs/pi/noonan/ak2267/LiverTF/fastq/fastqc
mkdir /gpfs/gibbs/pi/noonan/ak2267/LiverTF/bam
mkdir /gpfs/gibbs/pi/noonan/ak2267/LiverTF/bam/log
mkdir /gpfs/gibbs/pi/noonan/ak2267/LiverTF/peaks
mkdir /gpfs/gibbs/pi/noonan/ak2267/LiverTF/peaks/log
mkdir /gpfs/gibbs/pi/noonan/ak2267/LiverTF/visualize
mkdir /gpfs/gibbs/pi/noonan/ak2267/LiverTF/visualize/log

##################
# MAKE SQ FILES FOR DOWNLOAD, ALIGN, CALL PEAKS
cd /home/ak2267/project/LiverTF/batch

# get summary text files
# https://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-1509/
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-1509/E-MTAB-1509.sdrf.txt

python prepFiles/prepFiles_processLiverTF.py E-MTAB-1509.sdrf.txt
# makes
# 220720_downloadFastq_jobFile.txt
# 220720_fastqc_jobFile.txt
# 220720_align_jobFile.txt
# 220720_sort_jobFile.txt

# 220721_combine_jobFile.txt
# 220721_filter_jobFile.txt
# 220721_rename_jobFile.txt
# 220721_callPeaks_jobFile.txt
# 220721_visualize_jobFile.txt

# DOWNLOAD to /gpfs/gibbs/pi/noonan/ak2267/LiverTF/fastq
dsq --job-file 220720_downloadFastq_jobFile.txt --mem-per-cpu 5G -c 1 --mail-type FAIL,END
sbatch dsq-220720_downloadFastq_jobFile-2022-07-20.sh # 15657188

# FASTQC
dsq --job-file 220720_fastqc_jobFile.txt --mem-per-cpu 5G -c 1 --mail-type FAIL,END
sbatch dsq-220720_fastqc_jobFile-2022-07-20.sh # 15657242

# ALIGN to /gpfs/gibbs/pi/noonan/ak2267/LiverTF/bam
dsq --job-file 220720_align_jobFile.txt --mem-per-cpu 5G -c 8 --mail-type FAIL,END
sbatch dsq-220720_align_jobFile-2022-07-20.sh # 15657296

# DELETE FASTQ
rm /gpfs/gibbs/pi/noonan/ak2267/LiverTF/fastq/*.fastq.gz

# SORT
dsq --job-file 220720_sort_jobFile.txt --mem-per-cpu 5G -c 2 --mail-type FAIL,END
sbatch dsq-220720_sort_jobFile-2022-07-20.sh # 15657351
rm /gpfs/gibbs/pi/noonan/ak2267/LiverTF/bam/*_bowtie2.bam

# COMBINE
cd /home/ak2267/project/LiverTF/batch
dsq --job-file 220720_combine_jobFile.txt --mem-per-cpu 5G -c 1 --mail-type FAIL,END
sbatch dsq-220720_combine_jobFile-2022-07-20.sh # 15657405
rm /gpfs/gibbs/pi/noonan/ak2267/LiverTF/bam/*_a_sorted.ba*
rm /gpfs/gibbs/pi/noonan/ak2267/LiverTF/bam/*_b_sorted.ba*
rm /gpfs/gibbs/pi/noonan/ak2267/LiverTF/bam/*_c_sorted.ba*
rm /gpfs/gibbs/pi/noonan/ak2267/LiverTF/bam/*_d_sorted.ba*
rm /gpfs/gibbs/pi/noonan/ak2267/LiverTF/bam/*_e_sorted.ba*

# FILTER
dsq --job-file 220720_filter_jobFile.txt --mem-per-cpu 5G -c 2 --mail-type FAIL,END
sbatch dsq-220720_filter_jobFile-2022-07-20.sh # 15657442 STOPPED HERE 7/20/22
rm /gpfs/gibbs/pi/noonan/ak2267/LiverTF/bam/*_sorted.bam*
rm /gpfs/gibbs/pi/noonan/ak2267/LiverTF/bam/*_intermediateFiltered.bam*
rm /gpfs/gibbs/pi/noonan/ak2267/LiverTF/bam/*.bai

# just run the contents of 220721_rename_jobFile.txt on interactive node - simple mv commands
# only renames IPs - inputs remain with individual names for clarity on next steps
cd /gpfs/gibbs/pi/noonan/ak2267/LiverTF/bam/

# CALL PEAKS
cd /home/ak2267/project/LiverTF/batch
dsq --job-file 220721_callPeaks_jobFile.txt --mem-per-cpu 5G -c 1 --mail-type FAIL,END
sbatch dsq-220721_callPeaks_jobFile-2022-07-21.sh # 15658402

# VISUALIZE
dsq --job-file 220721_visualize_jobFile.txt --mem-per-cpu 5G -c 1 --mail-type FAIL,END
sbatch dsq-220721_visualize_jobFile-2022-07-21.sh # 15658880

# Look at files on browser
ssh ak2267@10.5.37.220
# cat new trackDb files onto old ones with Roller data and transfer to server
for species in rheMac10 mm39 rn7 canFam6
do
    cat /home/ak2267/project/Roller/ChIP/batch/trackDb_${species}.txt trackDb_${species}_LiverTFs.txt > trackDb_${species}.txt
    scp trackDb_${species}.txt ak2267@10.5.37.220:/home/ak2267/akocher_www/CGI/${species}
done

# transfer bw and bb files
cd /gpfs/gibbs/pi/noonan/ak2267/LiverTF/visualize/
for species in rheMac10 mm39 rn7 canFam6
do
    for mark in CEBPA FOXA1 HNF4A HNF6 inputTF
    do
	for i in $(ls ${species}*${mark}*)
	do
	    scp ${i} ak2267@10.5.37.220:/home/ak2267/akocher_www/CGI/${species}
	done
	sleep 10s
done;done

# RESTRICT TO REPLICATING PEAKS and write files into project space
# follow same decisions about which replicates to exclude as Ballester 2014 eLife, based on "substantially less peaks"
# Exclude lower replicate for: rheMac10 CEBPA rep 1, rheMac10 HNF4A, rheMac10 FOXA1, canFam6 ONECUT1
# These lower replicates aren't even available for download from ArrayExpress, so I only have one replicate for all of these

# I note that the following replicates also have very few peaks:
# canFam6 FOXA1 rep 1 --> remove before intersections
# rheMac10 HNF6 rep 1 --> remove before intersections
# rn7 HNF6 rep 2 --> remove before intersections
# mm39 FOXA1 rep 2 --> remove before intersections

# and these are borderline and may be removed from downstream analysis: rheMac10 CEBPA rep 1 (the only rep), canFam6 HNF6 rep 1 (the only rep)

mkdir /gpfs/gibbs/pi/noonan/ak2267/LiverTF/peaks/intersection
cd /gpfs/gibbs/pi/noonan/ak2267/LiverTF/peaks
for species in mm39 rn7
do
    bedtools intersect -a ${species}_CEBPA_1_peaks.narrowPeak -b ${species}_CEBPA_2_peaks.narrowPeak > intersection/${species}_liver_CEBPA_intersection.bed
    bedtools intersect -a ${species}_HNF4A_1_peaks.narrowPeak -b ${species}_HNF4A_2_peaks.narrowPeak > intersection/${species}_liver_HNF4A_intersection.bed
done

species=mm39
bedtools intersect -a ${species}_HNF6_1_peaks.narrowPeak -b ${species}_HNF6_2_peaks.narrowPeak > intersection/${species}_liver_HNF6_intersection.bed
cp ${species}_FOXA1_1_peaks.narrowPeak intersection/${species}_liver_FOXA1_intersection.bed

species=rn7
cp ${species}_HNF6_1_peaks.narrowPeak intersection/${species}_liver_HNF6_intersection.bed
bedtools intersect -a ${species}_FOXA1_1_peaks.narrowPeak -b ${species}_FOXA1_2_peaks.narrowPeak > intersection/${species}_liver_FOXA1_intersection.bed

species=rheMac10
cp ${species}_CEBPA_1_peaks.narrowPeak intersection/${species}_liver_CEBPA_intersection.bed
cp ${species}_FOXA1_1_peaks.narrowPeak intersection/${species}_liver_FOXA1_intersection.bed
cp ${species}_HNF4A_1_peaks.narrowPeak intersection/${species}_liver_HNF4A_intersection.bed
cp ${species}_HNF6_2_peaks.narrowPeak intersection/${species}_liver_HNF6_intersection.bed

species=canFam6
bedtools intersect -a ${species}_CEBPA_1_peaks.narrowPeak -b ${species}_CEBPA_2_peaks.narrowPeak > intersection/${species}_liver_CEBPA_intersection.bed
cp ${species}_FOXA1_2_peaks.narrowPeak intersection/${species}_liver_FOXA1_intersection.bed
bedtools intersect -a ${species}_HNF4A_1_peaks.narrowPeak -b ${species}_HNF4A_2_peaks.narrowPeak > intersection/${species}_liver_HNF4A_intersection.bed
cp ${species}_HNF6_1_peaks.narrowPeak intersection/${species}_liver_HNF6_intersection.bed

# Peaks are here:
/gpfs/gibbs/pi/noonan/ak2267/LiverTF/peaks/intersection/${species}_liver_${TF}_intersection.bed

# write to repsToUseLiverTF.txt
echo -e rheMac10'\t'liver'\t'CEBPA'\t'1 >> /gpfs/gibbs/pi/noonan/ak2267/LiverTF/peaks/repsToUseLiverTF.txt
echo -e rheMac10'\t'liver'\t'FOXA1'\t'1 >> /gpfs/gibbs/pi/noonan/ak2267/LiverTF/peaks/repsToUseLiverTF.txt
echo -e rheMac10'\t'liver'\t'HNF4A'\t'1 >> /gpfs/gibbs/pi/noonan/ak2267/LiverTF/peaks/repsToUseLiverTF.txt
echo -e rheMac10'\t'liver'\t'HNF6'\t'2 >> /gpfs/gibbs/pi/noonan/ak2267/LiverTF/peaks/repsToUseLiverTF.txt
for species in mm39 rn7
do
    for TF in CEBPA HNF4A
    do
	echo -e ${species}'\t'liver'\t'${TF}'\t'1,2 >> /gpfs/gibbs/pi/noonan/ak2267/LiverTF/peaks/repsToUseLiverTF.txt
    done;done >> /gpfs/gibbs/pi/noonan/ak2267/LiverTF/peaks/repsToUseLiverTF.txt
echo -e mm39'\t'liver'\t'HNF6'\t'1,2 >> /gpfs/gibbs/pi/noonan/ak2267/LiverTF/peaks/repsToUseLiverTF.txt
echo -e rn7'\t'liver'\t'HNF6'\t'1 >> /gpfs/gibbs/pi/noonan/ak2267/LiverTF/peaks/repsToUseLiverTF.txt
echo -e mm39'\t'liver'\t'FOXA1'\t'1 >> /gpfs/gibbs/pi/noonan/ak2267/LiverTF/peaks/repsToUseLiverTF.txt
echo -e rn7'\t'liver'\t'FOXA1'\t'1,2 >> /gpfs/gibbs/pi/noonan/ak2267/LiverTF/peaks/repsToUseLiverTF.txt

echo -e canFam6'\t'liver'\t'CEBPA'\t'1,2 >> /gpfs/gibbs/pi/noonan/ak2267/LiverTF/peaks/repsToUseLiverTF.txt
echo -e canFam6'\t'liver'\t'FOXA1'\t'2 >> /gpfs/gibbs/pi/noonan/ak2267/LiverTF/peaks/repsToUseLiverTF.txt
echo -e canFam6'\t'liver'\t'HNF4A'\t'1,2 >> /gpfs/gibbs/pi/noonan/ak2267/LiverTF/peaks/repsToUseLiverTF.txt
echo -e canFam6'\t'liver'\t'HNF6'\t'1 >> /gpfs/gibbs/pi/noonan/ak2267/LiverTF/peaks/repsToUseLiverTF.txt

