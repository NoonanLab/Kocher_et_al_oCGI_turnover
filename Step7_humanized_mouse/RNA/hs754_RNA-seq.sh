# 10/26/22
# RNA-seq analysis from hs754 humanized mouse

# Purpose: map, do differential expression analysis
# Samples: diencephalon from e11.5 & e17.5, 4 replicates per genotype (WT + HUM)

# work here:
mkdir /home/ak2267/project/hs754_RNA
cd /home/ak2267/project/hs754_RNA

mkdir /home/ak2267/project/hs754_RNA/fastqc
mkdir /home/ak2267/scratch60/hs754_RNA/
mkdir /home/ak2267/scratch60/hs754_RNA/bam

# sample info is in hs754_RNA-seq_sampleInfo.txt

# MAKE STAR INDEX FOR mm39 using Gencode vM31 on GRCm39
mkdir /home/ak2267/project/hs754_RNA/STARindex
cd /home/ak2267/project/hs754_RNA/STARindex
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M31/gencode.vM31.basic.annotation.gtf.gz

echo 'cd /home/ak2267/project/hs754_RNA/STARindex ; source /home/ak2267/.bashrc ; source /home/ak2267/.bash_profile ; module load STAR/2.7.9a-GCCcore-10.2.0 ; STAR --runThreadN 8 --runMode genomeGenerate --genomeDir /home/ak2267/project/hs754_RNA/STARindex/mm39 --genomeFastaFiles /gpfs/gibbs/pi/noonan/ak2267/genomes/mm39.fa --sjdbGTFfile /home/ak2267/project/hs754_RNA/STARindex/gencode.vM31.basic.annotation.gtf --sjdbOverhang 149' > 221027_makeSTARindex_jobFile.txt

dsq --job-file 221027_makeSTARindex_jobFile.txt --mem-per-cpu 5G -c 8 --mail-type FAIL,END
sbatch dsq-221027_makeSTARindex_jobFile-2022-10-27.sh # 18083447

# MAKE JOB FILES with /home/ak2267/project/hs754_RNA/prepFiles_RNA-seq.py
python prepFiles_hs754_RNA-seq.py hs754_RNA-seq_sampleInfo.txt

# makes the following job files:
# 221027_fastqc_jobFile.txt
# 221027_align_jobFile.txt

# run FASTQC
dsq --job-file 221027_fastqc_jobFile.txt --mem-per-cpu 5G -c 1 --mail-type FAIL,END
sbatch dsq-221027_fastqc_jobFile-2022-10-27.sh # 18084380

# run ALIGN which also does COUNT
dsq --job-file 221027_align_jobFile.txt --mem-per-cpu 5G -c 8 --mail-type FAIL,END
sbatch dsq-221027_align_jobFile-2022-10-27.sh # 18085231

# COUNT FILES are in format 
# F_e11.5_WT_1ReadsPerGene.out.tab
# where col 1 is ENSG and col4 is reverse strand counts (what we want to use)

# gzip and download
cd /home/ak2267/scratch60/hs754_RNA/bam/
mkdir results_hs754_RNA-seq
mv *ReadsPerGene.out.tab results_hs754_RNA-seq
tar -zcvf results_hs754_RNA-seq.gz results_hs754_RNA-seq/
cd /Users/acadiak/Desktop/CGI/hs754/RNA
scp ak2267@ruddle.hpc.yale.edu:/home/ak2267/scratch60/hs754_RNA/bam/results_hs754_RNA-seq.gz .

# analyze in R (hs754_RNA-seq.R)





