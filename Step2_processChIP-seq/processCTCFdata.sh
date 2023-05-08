# 6/15/22
# Purpose: download & process CTCF data from Schmidt et al 2012 Cell

# Data exists for CTCF in liver from 6 mammals: Human, Rhesus, Mouse, Rat, Dog, Opossum
# Note opossum data is a custom antibody different than the others (CTCF is apparently different in marsupials), and human is hepatocytes not whole liver
# I will use a subset that matches the species in Roller, Stamper et al 2021: Rhesus, Mouse, Rat, Dog

# work here
cd /home/ak2267/project/CTCF/batch

##################
# MAKE SQ FILES FOR DOWNLOAD, ALIGN, CALL PEAKS

# get summary text files
# https://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-437/
# https://www.ebi.ac.uk/ena/browser/view/PRJEB2329?show=reads
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-437/E-MTAB-437.sdrf.txt

python prepFiles_processCTCF.py E-MTAB-437.sdrf.txt
# makes
# 220620_downloadFastq_jobFile.txt
# 220620_fastqc_jobFile.txt
# 220620_align_jobFile.txt
# 220620_sort_jobFile.txt
# 220620_combine_jobFile.txt
# 220620_filter_jobFile.txt
# 220620_callPeaks_jobFile.txt
# 220724_visualize_jobFile.txt

# DOWNLOAD to /home/ak2267/scratch60/CTCF/fastq
dsq --job-file 220620_downloadFastq_jobFile.txt --mem-per-cpu 5G -c 1 --mail-type FAIL,END
sbatch dsq-220620_downloadFastq_jobFile-2022-06-20.sh # 15512783

# FASTQC
dsq --job-file 220620_fastqc_jobFile.txt --mem-per-cpu 5G -c 1 --mail-type FAIL,END
sbatch dsq-220620_fastqc_jobFile-2022-06-20.sh # 15512810

# ALIGN to /home/ak2267/scratch60/bam
dsq --job-file 220620_align_jobFile.txt --mem-per-cpu 5G -c 8 --mail-type FAIL,END
sbatch dsq-220620_align_jobFile-2022-06-20.sh # 15512898

# DELETE FASTQ
rm /home/ak2267/scratch60/CTCF/fastq/*.fastq.gz

# SORT
dsq --job-file 220620_sort_jobFile.txt --mem-per-cpu 5G -c 2 --mail-type FAIL,END
sbatch dsq-220620_sort_jobFile-2022-06-20.sh # 15512970
rm /home/ak2267/scratch60/CTCF/bam/*_bowtie2.bam

# COMBINE
dsq --job-file 220620_combine_jobFile.txt --mem-per-cpu 5G -c 1 --mail-type FAIL,END
sbatch dsq-220620_combine_jobFile-2022-06-20.sh # 15513065
rm /home/ak2267/scratch60/CTCF/bam/*_a_sorted.bam*
rm /home/ak2267/scratch60/CTCF/bam/*_b_sorted.bam*
rm /home/ak2267/scratch60/CTCF/bam/*_c_sorted.bam*

# FILTER
dsq --job-file 220620_filter_jobFile.txt --mem-per-cpu 5G -c 2 --mail-type FAIL,END
sbatch dsq-220620_filter_jobFile-2022-06-20.sh # 15513097
rm /home/ak2267/scratch60/CTCF/bam/*_sorted.bam*
rm /home/ak2267/scratch60/CTCF/bam/*_intermediateFiltered.bam*

# CALL PEAKS
dsq --job-file 220620_callPeaks_jobFile.txt --mem-per-cpu 5G -c 1 --mail-type FAIL,END
sbatch dsq-220620_callPeaks_jobFile-2022-06-20.sh # 15513178

# on 7/24/22 move files into the file structure in /gpfs/gibbs/pi/noonan/ak2267/LiverTF/ so they aren't deleted
mkdir /gpfs/gibbs/pi/noonan/ak2267/LiverTF/bam/CTCF
mv ~/scratch60/CTCF/bam/* /gpfs/gibbs/pi/noonan/ak2267/LiverTF/bam/CTCF
mkdir /gpfs/gibbs/pi/noonan/ak2267/LiverTF/peaks/CTCF
mv ~/scratch60/CTCF/peaks/* /gpfs/gibbs/pi/noonan/ak2267/LiverTF/peaks/CTCF

# VISUALIZE
cd /home/ak2267/project/CTCF/batch
python prepFiles/prepFiles_processCTCF.py E-MTAB-437.sdrf.txt
dsq --job-file 220724_visualize_jobFile.txt --mem-per-cpu 5G -c 2 --mail-type FAIL,END
sbatch dsq-220724_visualize_jobFile-2022-07-24.sh # 15675936

# move to browser - have to rename CTCF inputs so they don't conflict with Roller input names
cd /gpfs/gibbs/pi/noonan/ak2267/LiverTF/visualize
for species in rheMac10 mm39 rn7 canFam6
do
    mv ${species}_liver_input_1.bw ${species}_liver_inputCTCF_1.bw
    mv ${species}_liver_input_2.bw ${species}_liver_inputCTCF_2.bw
    scp ${species}_liver*CTCF* ak2267@10.5.37.220:/home/ak2267/akocher_www/CGI/${species}
done

### 7/27/22 
# Make replicating peaks

mkdir /gpfs/gibbs/pi/noonan/ak2267/LiverTF/peaks/CTCF/intersection
cd /gpfs/gibbs/pi/noonan/ak2267/LiverTF/peaks/CTCF

for species in rheMac10 mm39 rn7 canFam6
do
    bedtools intersect -a ${species}_liver_CTCF_1_peaks.narrowPeak -b ${species}_liver_CTCF_2_peaks.narrowPeak > intersection/${species}_liver_CTCF_intersection.bed
done
#cp canFam6_liver_CTCF_2_peaks.narrowPeak intersection/canFam6_liver_CTCF_intersection.bed
# intersection peaks are here
/gpfs/gibbs/pi/noonan/ak2267/LiverTF/peaks/CTCF/intersection/${species}_liver_CTCF_intersection.bed

# write to repsToUseLiverTF.txt
for species in rheMac10 mm39 rn7 canFam6
do
    echo -e $species'\t'liver'\t'CTCF'\t'1,2
done >> /gpfs/gibbs/pi/noonan/ak2267/LiverTF/peaks/repsToUseLiverTF.txt

#echo -e canFam6'\t'liver'\t'CTCF'\t'2 >> /gpfs/gibbs/pi/noonan/ak2267/LiverTF/peaks/repsToUseLiverTF.txt
