# 11/15/22
# Purpose: process Roller RNA and analyze effect of species-specific active CGIs on gene expression

# make STAR indexes for Roller genomes in UCSC coordinates
# have to redo from 220510_RollerRNA.sh because they were in scratch
# need to run STAR with these indexes, and have it output counts based on ENSEMBL GTFs but translated to UCSC chr names

# Work here
cd /home/ak2267/project/Roller/RNA/221121_RNA

# Download metadata file for Roller data
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-8122/E-MTAB-8122.sdrf.txt

# Download ENSEMBL GTFs (see Genome_and_annotation_info.xlsx for info on these files
# all are based on the same NCBI assemblies as the UCSC genomes I used for ChIP data
cd /home/ak2267/genomes/ENSEMBL
wget http://ftp.ensembl.org/pub/current_gtf/macaca_mulatta/Macaca_mulatta.Mmul_10.106.gtf.gz
wget http://ftp.ensembl.org/pub/current_gtf/mus_musculus/Mus_musculus.GRCm39.106.gtf.gz
wget http://ftp.ensembl.org/pub/current_gtf/rattus_norvegicus/Rattus_norvegicus.mRatBN7.2.106.gtf.gz
wget http://ftp.ensembl.org/pub/current_gtf/sus_scrofa/Sus_scrofa.Sscrofa11.1.106.gtf.gz
wget http://ftp.ensembl.org/pub/current_gtf/felis_catus/Felis_catus.Felis_catus_9.0.106.gtf.gz
wget http://ftp.ensembl.org/pub/current_gtf/equus_caballus/Equus_caballus.EquCab3.0.106.gtf.gz
gunzip *.gtf.gz

####### MAKE ENSEMBL GTF FILES BUT WITH CHROMOSOMES NAMED AS IN UCSC FASTA FILES

## Assess chromosome names in ENSEMBL vs UCSC files
# example for rhesus
cat Macaca_mulatta.Mmul_10.106.gtf | cut -f 1 | sort -u | wc -l # 334
cut -f 1 ~/genomes/RefSeq/UCSC_RefSeq/rheMac10.refGene.gtf | sort -u | wc -l #46
# UCSC collapses a lot of the unassembled scaffolds

## python script to convert ENSEMBL chromosome names to UCSC
# simply add 'chr' in front of chromosome numbers (if number is 1-2 digits) - otherwise remove from GTF file
cd /home/ak2267/project/Roller/RNA/221121_RNA
python convertChrNamesToUCSCinGTF.py /home/ak2267/genomes/ENSEMBL/Equus_caballus.EquCab3.0.106.gtf GTF/equCab3.ENSEMBLinUCSC.gtf
python convertChrNamesToUCSCinGTF.py /home/ak2267/genomes/ENSEMBL/Felis_catus.Felis_catus_9.0.106.gtf GTF/felCat9.ENSEMBLinUCSC.gtf
python convertChrNamesToUCSCinGTF.py /home/ak2267/genomes/ENSEMBL/Macaca_mulatta.Mmul_10.106.gtf GTF/rheMac10.ENSEMBLinUCSC.gtf
python convertChrNamesToUCSCinGTF.py /home/ak2267/genomes/ENSEMBL/Mus_musculus.GRCm39.106.gtf GTF/mm39.ENSEMBLinUCSC.gtf
python convertChrNamesToUCSCinGTF.py /home/ak2267/genomes/ENSEMBL/Rattus_norvegicus.mRatBN7.2.106.gtf GTF/rn7.ENSEMBLinUCSC.gtf
python convertChrNamesToUCSCinGTF.py /home/ak2267/genomes/ENSEMBL/Sus_scrofa.Sscrofa11.1.106.gtf GTF/susScr11.ENSEMBLinUCSC.gtf

####### OVERVIEW OF 1:1 ORTHOLOG FILES

# Download 1:1 ortholog files from BioMart into /Users/acadiak/Desktop/H3K4me3/Roller/RNA/1-1_Orthologs
# Filter --> Multi species comparisons --> check "Homologue filters" and select other genome.
# Then add other genome's Gene Stable ID under Attributes --> Homologues (Max select 6 orthologs)
# Download all files as .txt.gz and upload to Ruddle
cd /Users/acadiak/Desktop/H3K4me3/Roller/RNA/
scp -r 1-1_Orthologs ak2267@ruddle.hpc.yale.edu:/home/ak2267/genomes/ENSEMBL/
cd /home/ak2267/genomes/ENSEMBL/

# Count 1:1 orthologs that fall on unassembled contigs to assess the impact of excluding them
python assessOrthologDistribution.py 1-1_Orthologs/mouse_rat.txt Mus_musculus.GRCm39.106.gtf 1 # 22646/22686 orthologs are on main chromosomes
# similar in all files so I think it's ok to exclude those not on main chromosomes

# Count 1:1 orthologs that have only one ortholog in the other species
python assessNumberOfOrthologs.py 1-1_Orthologs/rhesus_mouse.txt # 15937/18448 orthologs are in fact 1:1
# similar in all files so I think it's ok to exclude those with 1:many or many:1 or many:many

####### MAKE REGULATORY DOMAIN BED FILES - EXTEND FROM TSSs IN ENSEMBL GTF FILES

# Assess biotypes & see if they are all captured by my other scripts for NCBI RefSeq files -> they are not
# Re-make assessment of biotypes
# see here for Biotype explanation https://www.gencodegenes.org/pages/biotypes.html
cd /home/ak2267/project/Roller/RNA/221121_RNA
for file in Equus_caballus.EquCab3.0.106.gtf Mus_musculus.GRCm39.106.gtf Felis_catus.Felis_catus_9.0.106.gtf Rattus_norvegicus.mRatBN7.2.106.gtf Macaca_mulatta.Mmul_10.106.gtf Sus_scrofa.Sscrofa11.1.106.gtf
do
    python outputUniqueBiotype.py /home/ak2267/genomes/ENSEMBL/${file}
done | cut -f 1 | sort -u
# see Genome_and_annotation_info.xlsx for summary (and RefSeq_summary.xlsx for details from doing this with NCBI & UCSC files)

##### Make proximal regulatory region bed files from translated GTFs (protein-coding genes only)
# Make files with 5kb upstream and 1kb downstream from every protein-coding gene
# Use a modified version of my previous script (makeCategoryBedFiles_plusIntrons.py)
# that additionally translates chromosome names (adding 'chr' to ENSEMBL names)
cd /home/ak2267/project/Roller/RNA/221121_RNA
python makeCategoryBedFiles_ENSEMBL.py /home/ak2267/genomes/ENSEMBL/Equus_caballus.EquCab3.0.106.gtf equCab3
python makeCategoryBedFiles_ENSEMBL.py /home/ak2267/genomes/ENSEMBL/Mus_musculus.GRCm39.106.gtf mm39
python makeCategoryBedFiles_ENSEMBL.py /home/ak2267/genomes/ENSEMBL/Felis_catus.Felis_catus_9.0.106.gtf felCat9
python makeCategoryBedFiles_ENSEMBL.py /home/ak2267/genomes/ENSEMBL/Rattus_norvegicus.mRatBN7.2.106.gtf rn7
python makeCategoryBedFiles_ENSEMBL.py /home/ak2267/genomes/ENSEMBL/Macaca_mulatta.Mmul_10.106.gtf rheMac10
python makeCategoryBedFiles_ENSEMBL.py /home/ak2267/genomes/ENSEMBL/Sus_scrofa.Sscrofa11.1.106.gtf susScr11

for species in rheMac10 mm39 rn7 susScr11 felCat9 equCab3
do
    echo $species
    sort -u ${species}_1_ENSEMBL_proteinCoding_promoters_preMerge.bed | sort -k1,1 -k2,2n > promoterBed/${species}_1_ENSEMBL_proteinCoding_promoters_sorted.bed
done
rm *_preMerge.bed


# MAKE STAR INDEXES
mkdir /home/ak2267/scratch60/Roller/genomes
mkdir /home/ak2267/scratch60/Roller/genomes/STARindexes
cd /home/ak2267/project/Roller/RNA/221121_RNA
for species in rheMac10 mm39 rn7 susScr11 felCat9 equCab3
do
    echo 'cd /home/ak2267/scratch60/Roller/genomes/STARindexes ; source /home/ak2267/.bashrc ; source /home/ak2267/.bash_profile ; module load STAR/2.7.9a-GCCcore-10.2.0 ; STAR --runThreadN 8 --runMode genomeGenerate --genomeDir '${species}' --genomeFastaFiles /gpfs/gibbs/pi/noonan/ak2267/genomes/'${species}'.fa --sjdbGTFfile /home/ak2267/project/Roller/RNA/221121_RNA/GTF/'${species}'.ENSEMBLinUCSC.gtf --sjdbOverhang 149'
done > 221122_makeSTARindexes_jobFile.txt

dsq --job-file 221122_makeSTARindexes_jobFile.txt --mem-per-cpu 5G -c 8 --mail-type FAIL,END
sbatch dsq-221122_makeSTARindexes_jobFile-2022-11-22.sh # 18394530

# one failed - rerun
dsqa -j 18394530 -s FAILED -f 221122_makeSTARindexes_jobFile.txt > 221122_makeSTARindexes_jobFile_redo1.txt
dsq --job-file 221122_makeSTARindexes_jobFile_redo1.txt --mem-per-cpu 5G -c 8 --mail-type FAIL,END
sbatch dsq-221122_makeSTARindexes_jobFile_redo1-2022-11-22.sh # 18394539

######## MAKE JOB FILES
cd /home/ak2267/project/Roller/RNA/221121_RNA

# Run prepFiles_RollerRNA_Nov22.py to make all job files for processing
python prepFiles_processRollerRNA_Nov22.py E-MTAB-8122.sdrf.txt
# makes
221122_downloadFastq_jobFile.txt
221122_fastqc_jobFile.txt
221122_align_jobFile.txt

# make folders for data
mkdir /home/ak2267/scratch60/Roller/fastq_RNA
mkdir /home/ak2267/scratch60/Roller/fastq_RNA/fastqc
mkdir /home/ak2267/scratch60/Roller/bam_RNA

# run DOWNLOAD
dsq --job-file 221122_downloadFastq_jobFile.txt --mem-per-cpu 5G -c 1 --mail-type FAIL,END
sbatch dsq-221122_downloadFastq_jobFile-2022-11-22.sh # 18394952
# accidentally left out horse - run prepFiles again
grep equCab3 221122_downloadFastq_jobFile.txt > 221122_downloadFastq_jobFile_plusEquCab3.txt 
dsq --job-file 221122_downloadFastq_jobFile_plusEquCab3.txt --mem-per-cpu 5G -c 1 --mail-type FAIL,END
sbatch dsq-221122_downloadFastq_jobFile_plusEquCab3-2022-11-22.sh # 18398801

# run FASTQC
dsq --job-file 221122_fastqc_jobFile.txt --mem-per-cpu 5G -c 1 --mail-type FAIL,END
sbatch dsq-221122_fastqc_jobFile-2022-11-27.sh # 18424039

# run Multiqc to visualize all samples
module load MultiQC/1.10.1-foss-2020b-Python-3.8.6
multiqc /home/ak2267/scratch60/Roller
/Users/acadiak/Desktop/CGI/RNA

# run ALIGN (also does COUNT in same step)
dsq --job-file 221122_align_jobFile.txt --mem-per-cpu 5G -c 8 --mail-type FAIL,END
sbatch dsq-221122_align_jobFile-2022-11-22.sh # 18401920 



##### Connect oCGIs to their nearest genes on either side (within 1 Mb)
# get oCGI info from here because it includes final filtering steps:
# /gpfs/gibbs/pi/noonan/ak2267/speciesPairs/CGIcentric/Roller_summaryFiles/rheMac10_mm39_brain_H3K4me3.txt 

# make bed file with oCGI coordinates in both species - restrict to A-only CGIs with A-only peaks and B-only CGIs with B-only peaks
# this happens separately for each species pair x tissue x mark
# INPUTS:
# 1) /gpfs/gibbs/pi/noonan/ak2267/speciesPairs/CGIcentric/Roller_summaryFiles/rheMac10_mm39_brain_H3K4me3.txt 
# 2) /gpfs/gibbs/pi/noonan/ak2267/speciesPairs/consensusCGIs/rheMac10_mm39/rheMac10_mm39_speciesA_reconciledCGI.bed
# 3) /gpfs/gibbs/pi/noonan/ak2267/speciesPairs/consensusCGIs/rheMac10_mm39/rheMac10_mm39_speciesB_reconciledCGI.bed
# OUTPUTS:
# 1) bed file with all A-only CGIs with A-only peaks and B-only CGIs with B-only peaks, in A coordinates
# 1) bed file with all A-only CGIs with A-only peaks and B-only CGIs with B-only peaks, in B coordinates

for speciesPair in rheMac10_mm39 rheMac10_rn7 rheMac10_susScr11 rheMac10_felCat9 rheMac10_equCab3 mm39_rn7 mm39_susScr11 mm39_felCat9 mm39_equCab3 rn7_susScr11 rn7_felCat9 rn7_equCab3 susScr11_felCat9 susScr11_equCab3 felCat9_equCab3
do
    for tissue in brain liver muscle testis
    do
        for mark in H3K4me3 H3K27ac H3K4me1
        do
            echo $speciesPair $tissue $mark
            python outputCGIsAsBed.py /gpfs/gibbs/pi/noonan/ak2267/speciesPairs/CGIcentric/Roller_summaryFiles/${speciesPair}_${tissue}_${mark}.txt /gpfs/gibbs/pi/noonan/ak2267/speciesPairs/consensusCGIs/${speciesPair}/${speciesPair}_speciesA_reconciledCGI.bed /gpfs/gibbs/pi/noonan/ak2267/speciesPairs/consensusCGIs/${speciesPair}/${speciesPair}_speciesB_reconciledCGI.bed
        done; done; done
        
# sort
for speciesPair in rheMac10_mm39 rheMac10_rn7 rheMac10_susScr11 rheMac10_felCat9 rheMac10_equCab3 mm39_rn7 mm39_susScr11 mm39_felCat9 mm39_equCab3 rn7_susScr11 rn7_felCat9 rn7_equCab3 susScr11_felCat9 susScr11_equCab3 felCat9_equCab3
do
    for tissue in brain liver muscle testis
    do
        for mark in H3K4me3 H3K27ac H3K4me1
        do
            sort -k1,1 -k2,2n nearestGene/${speciesPair}_${tissue}_${mark}_coordInA.bed > nearestGene/${speciesPair}_${tissue}_${mark}_coordInA_sorted.bed
            sort -k1,1 -k2,2n nearestGene/${speciesPair}_${tissue}_${mark}_coordInB.bed > nearestGene/${speciesPair}_${tissue}_${mark}_coordInB_sorted.bed
        done; done; done
rm nearestGene/*A.bed
rm nearestGene/*B.bed

# connect to nearest genes within 1 Mb
for speciesPair in rheMac10_mm39 rheMac10_rn7 rheMac10_susScr11 rheMac10_felCat9 rheMac10_equCab3 mm39_rn7 mm39_susScr11 mm39_felCat9 mm39_equCab3 rn7_susScr11 rn7_felCat9 rn7_equCab3 susScr11_felCat9 susScr11_equCab3 felCat9_equCab3
do
    speciesA="$(echo $speciesPair | cut -d'_' -f 1)"
    speciesB="$(echo $speciesPair | cut -d'_' -f 2)"
    for tissue in brain liver muscle testis
    do
        for mark in H3K4me3 H3K27ac H3K4me1
        do
            echo $speciesPair $tissue $mark
            bedtools closest -id -D ref -a nearestGene/${speciesPair}_${tissue}_${mark}_coordInA_sorted.bed -b promoterBed/${speciesA}_1_ENSEMBL_proteinCoding_promoters_sorted.bed > nearestGene/${speciesPair}_${tissue}_${mark}_speciesA_nearestUpstream.txt
            bedtools closest -iu -D ref -a nearestGene/${speciesPair}_${tissue}_${mark}_coordInA_sorted.bed -b promoterBed/${speciesA}_1_ENSEMBL_proteinCoding_promoters_sorted.bed > nearestGene/${speciesPair}_${tissue}_${mark}_speciesA_nearestDownstream.txt
            bedtools closest -id -D ref -a nearestGene/${speciesPair}_${tissue}_${mark}_coordInB_sorted.bed -b promoterBed/${speciesB}_1_ENSEMBL_proteinCoding_promoters_sorted.bed > nearestGene/${speciesPair}_${tissue}_${mark}_speciesB_nearestUpstream.txt
            bedtools closest -iu -D ref -a nearestGene/${speciesPair}_${tissue}_${mark}_coordInB_sorted.bed -b promoterBed/${speciesB}_1_ENSEMBL_proteinCoding_promoters_sorted.bed > nearestGene/${speciesPair}_${tissue}_${mark}_speciesB_nearestDownstream.txt
    done; done; done

# rename 1:1 ortholog files to make them easier to reference in the unix loop below
cp /home/ak2267/genomes/ENSEMBL/1-1_Orthologs/rhesus_mouse.txt /home/ak2267/genomes/ENSEMBL/1-1_Orthologs/rheMac10_mm39.txt 
cp /home/ak2267/genomes/ENSEMBL/1-1_Orthologs/rhesus_rat.txt /home/ak2267/genomes/ENSEMBL/1-1_Orthologs/rheMac10_rn7.txt 
cp /home/ak2267/genomes/ENSEMBL/1-1_Orthologs/rhesus_pig.txt /home/ak2267/genomes/ENSEMBL/1-1_Orthologs/rheMac10_susScr11.txt 
cp /home/ak2267/genomes/ENSEMBL/1-1_Orthologs/rhesus_cat.txt /home/ak2267/genomes/ENSEMBL/1-1_Orthologs/rheMac10_felCat9.txt 
cp /home/ak2267/genomes/ENSEMBL/1-1_Orthologs/rhesus_horse.txt /home/ak2267/genomes/ENSEMBL/1-1_Orthologs/rheMac10_equCab3.txt 
cp /home/ak2267/genomes/ENSEMBL/1-1_Orthologs/mouse_rat.txt /home/ak2267/genomes/ENSEMBL/1-1_Orthologs/mm39_rn7.txt 
cp /home/ak2267/genomes/ENSEMBL/1-1_Orthologs/mouse_pig.txt /home/ak2267/genomes/ENSEMBL/1-1_Orthologs/mm39_susScr11.txt 
cp /home/ak2267/genomes/ENSEMBL/1-1_Orthologs/mouse_cat.txt /home/ak2267/genomes/ENSEMBL/1-1_Orthologs/mm39_felCat9.txt 
cp /home/ak2267/genomes/ENSEMBL/1-1_Orthologs/mouse_horse.txt /home/ak2267/genomes/ENSEMBL/1-1_Orthologs/mm39_equCab3.txt 
cp /home/ak2267/genomes/ENSEMBL/1-1_Orthologs/rat_pig.txt /home/ak2267/genomes/ENSEMBL/1-1_Orthologs/rn7_susScr11.txt 
cp /home/ak2267/genomes/ENSEMBL/1-1_Orthologs/rat_cat.txt /home/ak2267/genomes/ENSEMBL/1-1_Orthologs/rn7_felCat9.txt 
cp /home/ak2267/genomes/ENSEMBL/1-1_Orthologs/rat_horse.txt /home/ak2267/genomes/ENSEMBL/1-1_Orthologs/rn7_equCab3.txt 
cp /home/ak2267/genomes/ENSEMBL/1-1_Orthologs/pig_cat.txt /home/ak2267/genomes/ENSEMBL/1-1_Orthologs/susScr11_felCat9.txt 
cp /home/ak2267/genomes/ENSEMBL/1-1_Orthologs/pig_horse.txt /home/ak2267/genomes/ENSEMBL/1-1_Orthologs/susScr11_equCab3.txt 
cp /home/ak2267/genomes/ENSEMBL/1-1_Orthologs/cat_horse.txt /home/ak2267/genomes/ENSEMBL/1-1_Orthologs/felCat9_equCab3.txt 


######## UPDATE PIPELINE to include assessment of 1:1 orthology in human AND number of total enhancers in the region
cd /home/ak2267/project/Roller/RNA/221121_RNA


##### Human orthology
# Download 1:1 orthology files between each species and human
# from BioMart (Ensembl Genes 108) into /Users/acadiak/Desktop/CGI/RNA/1-1_Orthologs

# Filter --> Multi species comparisons --> check "Homologue filters" and select other genome.
# Then add other genome's Gene Stable ID under Attributes --> Homologues (Max select 6 orthologs)
# Download all files as .txt.gz
# Also download SINGLE file with ALL human genes and orthologs in the 6 species - try using this one, upload to Ruddle

cd /Users/acadiak/Desktop/CGI/RNA/
scp -r 1-1_Orthologs ak2267@ruddle.hpc.yale.edu:/home/ak2267/project/Roller/RNA/221121_RNA
cd /home/ak2267/project/Roller/RNA/221121_RNA/1-1_Orthologs
gunzip *.gz
mv *.txt /home/ak2267/genomes/ENSEMBL/1-1_Orthologs

# update script to incorporate human ortholog info
cd /home/ak2267/project/Roller/RNA/221121_RNA

##### Other peaks in region
mkdir allPeaks

for species in rheMac10 mm39 rn7 susScr11 felCat9 equCab3
do
    for tissue in brain liver muscle testis
    do
        for mark in H3K4me3 H3K27ac H3K4me1
        do
            echo $species $tissue $mark
            sort -k1,1 -k2,2n /gpfs/gibbs/pi/noonan/ak2267/Roller/peaks/intersection/${species}_${tissue}_${mark}_intersection.bed | bedtools intersect -v -a - -b /home/ak2267/genomes/RefSeq/featureAnnotations/${species}_allFeatures.bed | bedtools closest -id -D ref -a - -b promoterBed/${species}_1_ENSEMBL_proteinCoding_promoters_sorted.bed > allPeaks/${species}_${tissue}_${mark}_nearestUpstream.txt
            sort -k1,1 -k2,2n /gpfs/gibbs/pi/noonan/ak2267/Roller/peaks/intersection/${species}_${tissue}_${mark}_intersection.bed | bedtools intersect -v -a - -b /home/ak2267/genomes/RefSeq/featureAnnotations/${species}_allFeatures.bed | bedtools closest -iu -D ref -a - -b promoterBed/${species}_1_ENSEMBL_proteinCoding_promoters_sorted.bed > allPeaks/${species}_${tissue}_${mark}_nearestDownstream.txt
        done; done; done

# union of all other peaks tied to nearest gene
rm allPeaks/*_cat.bed
for species in rheMac10 mm39 rn7 susScr11 felCat9 equCab3
do
    for tissue in brain liver muscle testis
    do
        for mark in H3K4me3 H3K27ac H3K4me1
        do
            echo $species $tissue $mark
            sort -k1,1 -k2,2n /gpfs/gibbs/pi/noonan/ak2267/Roller/peaks/intersection/${species}_${tissue}_${mark}_intersection.bed | bedtools intersect -v -a - -b /home/ak2267/genomes/RefSeq/featureAnnotations/${species}_allFeatures.bed >> allPeaks/${species}_${tissue}_cat.bed
        done
        
        sort -k1,1 -k2,2n allPeaks/${species}_${tissue}_cat.bed > allPeaks/${species}_${tissue}_sort.bed
        bedtools merge -c 4 -o collapse -i allPeaks/${species}_${tissue}_sort.bed > allPeaks/${species}_${tissue}_merged.bed
                        
        bedtools closest -id -D ref -a allPeaks/${species}_${tissue}_merged.bed -b promoterBed/${species}_1_ENSEMBL_proteinCoding_promoters_sorted.bed > allPeaks/${species}_${tissue}_allMarks_nearestUpstream.txt
        bedtools closest -iu -D ref -a allPeaks/${species}_${tissue}_merged.bed -b promoterBed/${species}_1_ENSEMBL_proteinCoding_promoters_sorted.bed > allPeaks/${species}_${tissue}_allMarks_nearestDownstream.txt
        done; done        

###### UPDATE: Make regulatory domain bed files instead of approach using bedtools closest

# FIRST (done above)
# Make proximal regulatory region bed files from translated GTFs (protein-coding genes only)
# Make files with 5kb upstream and 1kb downstream from every protein-coding gene
# Use a modified version of my previous script (makeCategoryBedFiles_plusIntrons.py)
# that additionally translates chromosome names (adding 'chr' to ENSEMBL names)
cd /home/ak2267/project/Roller/RNA/221121_RNA
python makeCategoryBedFiles_ENSEMBL.py /home/ak2267/genomes/ENSEMBL/Equus_caballus.EquCab3.0.106.gtf equCab3
python makeCategoryBedFiles_ENSEMBL.py /home/ak2267/genomes/ENSEMBL/Mus_musculus.GRCm39.106.gtf mm39
python makeCategoryBedFiles_ENSEMBL.py /home/ak2267/genomes/ENSEMBL/Felis_catus.Felis_catus_9.0.106.gtf felCat9
python makeCategoryBedFiles_ENSEMBL.py /home/ak2267/genomes/ENSEMBL/Rattus_norvegicus.mRatBN7.2.106.gtf rn7
python makeCategoryBedFiles_ENSEMBL.py /home/ak2267/genomes/ENSEMBL/Macaca_mulatta.Mmul_10.106.gtf rheMac10
python makeCategoryBedFiles_ENSEMBL.py /home/ak2267/genomes/ENSEMBL/Sus_scrofa.Sscrofa11.1.106.gtf susScr11

for species in rheMac10 mm39 rn7 susScr11 felCat9 equCab3
do
    echo $species
    sort -u ${species}_1_ENSEMBL_proteinCoding_promoters_preMerge.bed | sort -k1,1 -k2,2n > promoterBed/${species}_1_ENSEMBL_proteinCoding_promoters_sorted.bed
done
rm *_preMerge.bed

# SECOND (done here)
# Extend proximal regulatory domains to nearest neighboring site, or up to 1 Mb
cd /home/ak2267/project/Roller/RNA/221121_RNA
mkdir regulatoryDomains
for species in rheMac10 mm39 rn7 susScr11 felCat9 equCab3
do
    echo $species
    python extendRegulatoryDomains.py promoterBed/${species}_1_ENSEMBL_proteinCoding_promoters_sorted.bed /home/ak2267/genomes/chrom.sizes/${species}.chrom.sizes > regulatoryDomains/${species}_regulatoryDomains_ENSEMBL.bed
done

###### Overlap CGIs with THESE regulatory domain files
mkdir overlappingRegulatoryDomains

for speciesPair in rheMac10_mm39 rheMac10_rn7 rheMac10_susScr11 rheMac10_felCat9 rheMac10_equCab3 mm39_rn7 mm39_susScr11 mm39_felCat9 mm39_equCab3 rn7_susScr11 rn7_felCat9 rn7_equCab3 susScr11_felCat9 susScr11_equCab3 felCat9_equCab3
do
    speciesA="$(echo $speciesPair | cut -d'_' -f 1)"
    speciesB="$(echo $speciesPair | cut -d'_' -f 2)"
    for tissue in brain liver muscle testis
    do
        for mark in H3K4me3 H3K27ac H3K4me1
        do
            echo $speciesPair $tissue $mark
            bedtools intersect -wao -a nearestGene/${speciesPair}_${tissue}_${mark}_coordInA_sorted.bed -b regulatoryDomains/${speciesA}_regulatoryDomains_ENSEMBL.bed > overlappingRegulatoryDomains/${speciesPair}_${tissue}_${mark}_speciesA.txt
            bedtools intersect -wao -a nearestGene/${speciesPair}_${tissue}_${mark}_coordInB_sorted.bed -b regulatoryDomains/${speciesB}_regulatoryDomains_ENSEMBL.bed > overlappingRegulatoryDomains/${speciesPair}_${tissue}_${mark}_speciesB.txt
    done; done; done
    
###### Overlap all peaks with THESE regulatory domain files

# individual peak types
for species in rheMac10 mm39 rn7 susScr11 felCat9 equCab3
do
    for tissue in brain liver muscle testis
    do
        for mark in H3K4me3 H3K27ac H3K4me1
        do
            echo $species $tissue $mark
            sort -k1,1 -k2,2n /gpfs/gibbs/pi/noonan/ak2267/Roller/peaks/intersection/${species}_${tissue}_${mark}_intersection.bed | bedtools intersect -v -a - -b /home/ak2267/genomes/RefSeq/featureAnnotations/${species}_allFeatures.bed | bedtools intersect -wao -a - -b regulatoryDomains/${species}_regulatoryDomains_ENSEMBL.bed > allPeaks/${species}_${tissue}_${mark}_intersectRegulatoryDomains.txt
        done; done; done

# union of all peaks tied to nearest gene
rm allPeaks/*_cat.bed
for species in rheMac10 mm39 rn7 susScr11 felCat9 equCab3
do
    for tissue in brain liver muscle testis
    do
        for mark in H3K4me3 H3K27ac H3K4me1
        do
            echo $species $tissue $mark
            sort -k1,1 -k2,2n /gpfs/gibbs/pi/noonan/ak2267/Roller/peaks/intersection/${species}_${tissue}_${mark}_intersection.bed | bedtools intersect -v -a - -b /home/ak2267/genomes/RefSeq/featureAnnotations/${species}_allFeatures.bed >> allPeaks/${species}_${tissue}_cat.bed
        done
        
        sort -k1,1 -k2,2n allPeaks/${species}_${tissue}_cat.bed > allPeaks/${species}_${tissue}_sort.bed
        bedtools merge -c 4 -o collapse -i allPeaks/${species}_${tissue}_sort.bed > allPeaks/${species}_${tissue}_merged.bed
                        
        bedtools intersect -wao -a allPeaks/${species}_${tissue}_merged.bed -b regulatoryDomains/${species}_regulatoryDomains_ENSEMBL.bed > allPeaks/${species}_${tissue}_intersectRegulatoryDomains.txt
        done; done        


########################################################################
##### Run script to make summary table by incorporating all of the above
for speciesPair in rheMac10_mm39 rheMac10_rn7 rheMac10_susScr11 rheMac10_felCat9 rheMac10_equCab3 mm39_rn7 mm39_susScr11 mm39_felCat9 mm39_equCab3 rn7_susScr11 rn7_felCat9 rn7_equCab3 susScr11_felCat9 susScr11_equCab3 felCat9_equCab3
do
    for tissue in brain liver muscle testis
    do
        for mark in H3K4me3 H3K27ac H3K4me1
        do
            echo $speciesPair $tissue $mark
            python makeTwoWayOrthologTable.py ${speciesPair} ${tissue} ${mark} ReadsPerGene.out.tab > ENSEMBL_orthologs/${speciesPair}_${tissue}_${mark}_ENSEMBLorthologs.txt
            done; done; done
            

###### gzip and download files
cd /home/ak2267/project/Roller/RNA/221121_RNA
tar -zcvf ENSEMBL_orthologs.gz ENSEMBL_orthologs/
cd /Users/acadiak/Desktop/CGI/RNA
scp ak2267@ruddle.hpc.yale.edu:/home/ak2267/project/Roller/RNA/221121_RNA/ENSEMBL_orthologs.gz .

            
# upload exon.lengths file (made in Fig6_RNA.R) to Ruddle
scp /Users/acadiak/Desktop/CGI/RNA/exon.lengths.txt ak2267@ruddle.hpc.yale.edu:/home/ak2267/project/Roller/RNA/221121_RNA/permutation
            
########################################################################
##### Run permutation test using the summary tables

for speciesPair in rheMac10_mm39 rheMac10_rn7 rheMac10_susScr11 rheMac10_felCat9 rheMac10_equCab3 mm39_rn7 mm39_susScr11 mm39_felCat9 mm39_equCab3 rn7_susScr11 rn7_felCat9 rn7_equCab3 susScr11_felCat9 susScr11_equCab3 felCat9_equCab3
do
    for tissue in brain liver muscle testis
    do
        for mark in H3K4me3 H3K27ac H3K4me1
        do
            echo 'cd /home/ak2267/project/Roller/RNA/221121_RNA/permutation ; source ~/.bashrc ; source ~/.bash_profile ; module load R ; Rscript RNA_resamplingTest.R '${speciesPair}'_'${tissue}'_'${mark}'_ENSEMBLorthologs.txt'
            done; done; done >> 230212_RNA_resamplingTest_jobFile.txt
            
dsq --job-file 230212_RNA_resamplingTest_jobFile.txt --mem-per-cpu 5G -c 1 --mail-type FAIL,END --max-jobs 60
sbatch dsq-230212_RNA_resamplingTest_jobFile-2023-02-12.sh # 21100275

# Summarize results of resampling in tables for use in R

for speciesPair in rheMac10_mm39 rheMac10_rn7 rheMac10_susScr11 rheMac10_felCat9 rheMac10_equCab3 mm39_rn7 mm39_susScr11 mm39_felCat9 mm39_equCab3 rn7_susScr11 rn7_felCat9 rn7_equCab3 susScr11_felCat9 susScr11_equCab3 felCat9_equCab3
do
    for tissue in brain liver muscle testis
    do
        for mark in H3K4me3 H3K27ac H3K4me1
        do
            python summarizeRNAresampling.py ${speciesPair} ${tissue} ${mark}
        done; done; done >> RNA_resamplingSummary.txt

# download this file
cd /Users/acadiak/Desktop/CGI/RNA
scp ak2267@ruddle.hpc.yale.edu:/home/ak2267/project/Roller/RNA/221121_RNA/permutation/RNA_resamplingSummary.txt .


# download file for example
cd /Users/acadiak/Desktop/CGI/RNA
scp ak2267@ruddle.hpc.yale.edu:/home/ak2267/project/Roller/RNA/221121_RNA/permutation/rn7_susScr11_brain_H3K27ac_A_resamplingMedians.txt .
scp ak2267@ruddle.hpc.yale.edu:/home/ak2267/project/Roller/RNA/221121_RNA/permutation/rn7_susScr11_brain_H3K27ac_B_resamplingMedians.txt .
scp ak2267@ruddle.hpc.yale.edu:/home/ak2267/project/Roller/RNA/221121_RNA/permutation/rn7_susScr11_brain_H3K27ac_observedMedians.txt .

