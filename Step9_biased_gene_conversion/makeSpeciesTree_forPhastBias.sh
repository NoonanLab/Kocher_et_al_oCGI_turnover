# 9/13/22
# Purpose: generate neutral.mod file for running phastBias
# using msa_view and phyloFit

# work here:
cd /gpfs/gibbs/pi/noonan/ak2267/speciesPairs/SequenceAnalysis/gBGC_analysis

# maf files (with my species extracted) are here:
/home/su57/scratch60/200120_120mammals/maf/chr*_subset.maf

##### STEP 1:
# run msa_view to extract 4-fold degenerate sites
# --4d: (For use with --features; assumes coding regions have feature type 'CDS')
#        Extract sufficient statistics for fourfold degenerate synonymous sites.
#        Implies --out-format SS --unordered-stats --tuple-size 3 --reverse-groups transcript_id.
# --features, -g <gff_fname>
#       (Requires --catmap) Read sequence annotations from the
#       specified file (GFF) and label the columns of the alignment
#       accordingly.  Note: UCSC BED and genepred formats are now
#       recognized as well.
# I have no idea what catmap means so try running without

# Download human GFF file from Gencode
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_41/gencode.v41.annotation.gff3.gz
gunzip gencode.v41.annotation.gff3.gz

# Split by chromosome and rename to match maf coordinates (chr1 -> hg38.chr1)
for num in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X
do
    echo $num
    head -6 gencode.v41.annotation.gff3 > gencode.v41.annotation_chr${num}.gff3
    grep -w chr${num} gencode.v41.annotation.gff3 | sed -e 's/chr'${num}'/hg38.chr'${num}'/g' >> gencode.v41.annotation_chr${num}.gff3
done
mkdir gencode; mv gencode.* gencode/

# Run msa_view
for num in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X
do
    echo 'cd /gpfs/gibbs/pi/noonan/ak2267/speciesPairs/SequenceAnalysis/gBGC_analysis ; source /home/ak2267/.bashrc ; source /home/ak2267/.bash_profile ; msa_view --order hg38,rheMac8,calJac3,mm10,rn6,susScr11,canFam3,felCat8,HLequCab3 --features gencode/gencode.v41.annotation_chr'${num}'.gff3 --4d /home/su57/scratch60/200120_120mammals/maf/chr'${num}'_subset.maf > chr'${num}'.ss'
done >> 220913_runMsaView_jobFile.txt
dsq --job-file 220913_runMsaView_jobFile.txt --mem-per-cpu 50G -c 1 --mail-type FAIL,END
sbatch dsq-220913_runMsaView_jobFile-2022-09-13.sh # 16707282

for num in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X
do
    mv hg38.chr${num}.ss chr${num}.ss
done

# Aggregate sufficient statistics across all chromosomes
msa_view --unordered-ss --out-format SS --tuple-size 3 --aggregate hg38,rheMac8,calJac3,mm10,rn6,susScr11,canFam3,felCat8,HLequCab3 chr1.ss chr2.ss chr3.ss chr4.ss chr5.ss chr6.ss chr7.ss chr8.ss chr9.ss chr10.ss chr11.ss chr12.ss chr13.ss chr14.ss chr15.ss chr16.ss chr17.ss chr18.ss chr19.ss chr20.ss chr21.ss chr22.ss chrX.ss > allChr.ss
msa_view allChr.ss --in-format SS --out-format SS --tuple-size 1 > allChr_tuple1.ss

##### STEP 2:
# run phyloFit 

phyloFit --tree "((((hg38,rheMac8),calJac3),(mm10,rn6)),((susScr11,(canFam3,felCat8)),HLequCab3))" --subst-mod REV --out-root AK_tree --msa-format SS allChr_tuple1.ss


##### STEP 3:
# run phastBias - done by SU 9/13/22
