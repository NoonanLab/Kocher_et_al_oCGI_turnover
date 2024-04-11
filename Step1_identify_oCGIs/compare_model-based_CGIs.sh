# 03/02/2024
# Purpose: compare the CGIs I used in the paper to the model-based CGIs defined
# by the Wu lab https://www.haowulab.org/software/makeCGI/index.html

# RESULTS ARE REVISED FIG S2

# Download model-based CGIs for species in my dataset
mkdir /home/ak2267/genomes/CGI/Model_based
cd /home/ak2267/genomes/CGI/Model_based

wget https://www.haowulab.org/software/makeCGI/model-based-cpg-islands-hg38.txt
wget https://www.haowulab.org/software/makeCGI/model-based-cpg-islands-hg19.txt
wget https://www.haowulab.org/software/makeCGI/model-based-cpg-islands-mm10.txt
wget https://www.haowulab.org/software/makeCGI/model-based-cpg-islands-mm9.txt
wget https://www.haowulab.org/software/makeCGI/model-based-cpg-islands-rn4.txt
wget https://www.haowulab.org/software/makeCGI/model-based-cpg-islands-rheMac2.txt
wget https://www.haowulab.org/software/makeCGI/model-based-cpg-islands-canFam2.txt
wget https://www.haowulab.org/software/makeCGI/model-based-cpg-islands-equCab2.txt

## NOTE THAT THESE FILES ARE 1-BASED AND NOT 0-BASED ^^^^ I fixed this on 3/5/2024


# My strategy for this analysis:
# Restrict both sets to oCGIs that meet the criteria in my study (intronic-intergenic), then compare by means of:
#   Simple overlaps to show how many WITHIN a genome overlap (stacked bars)

mkdir /gpfs/gibbs/pi/noonan/ak2267/Revisions/Model-based_CGIs
mkdir /gpfs/gibbs/pi/noonan/ak2267/Revisions/Model-based_CGIs/Simple_overlaps

cd /gpfs/gibbs/pi/noonan/ak2267/Revisions/Model-based_CGIs/Simple_overlaps
mkdir Bed



# found errors in rn4 file
# line 13702 contained a coordinate of "9.5e+07" instead of 95000000
# line 49394 contained a coordinate of "1.2e+08" instead of 120000000
# corrected using the following:
cd /home/ak2267/genomes/CGI/Model_based

head -13701 model-based-cpg-islands-rn4.txt > rn4_pt1.bed
head -13702 model-based-cpg-islands-rn4.txt | tail -1 > rn4_pt2.bed # manually edit in vim
tail -78098 model-based-cpg-islands-rn4.txt > rn4_pt3.bed
cat rn4_pt1.bed rn4_pt2.bed rn4_pt3.bed > rn4_cat.bed
mv rn4_cat.bed model-based-cpg-islands-rn4.txt

head -49393 model-based-cpg-islands-rn4.txt > rn4_pt1.bed
head -49394 model-based-cpg-islands-rn4.txt | tail -1 > rn4_pt2.bed # manually edit in vim
tail -42406 model-based-cpg-islands-rn4.txt > rn4_pt3.bed
cat rn4_pt1.bed rn4_pt2.bed rn4_pt3.bed > rn4_cat.bed
mv rn4_cat.bed model-based-cpg-islands-rn4.txt


# get chains I didn't have
cd /home/ak2267/genomes/chain/

wget https://hgdownload.soe.ucsc.edu/goldenPath/rn4/liftOver/rn4ToRn6.over.chain.gz
wget https://hgdownload.soe.ucsc.edu/goldenPath/rn6/liftOver/rn6ToRn7.over.chain.gz
wget https://hgdownload.soe.ucsc.edu/goldenPath/canFam2/liftOver/canFam2ToCanFam6.over.chain.gz
wget https://hgdownload.soe.ucsc.edu/goldenPath/equCab2/liftOver/equCab2ToEquCab3.over.chain.gz


#############################################################################################
###### GOAL: compare numbers within a species                                          ###### 
###### Now I am re-doing this with new files that are NOT merged within 200 bp         ###### 
###### And that have info on length, CpG numbers, CpG obs/exp in the name column       ###### 
###### I'm also including UCSC CGIs to compare them to my lists                        ###### 
#############################################################################################


##### STEP 0:
# DOWNLOAD UCSC CGI LISTS (distinct from UCSC_AL lists!)

mkdir /home/ak2267/genomes/CGI/Standard_UCSC
cd /home/ak2267/genomes/CGI/Standard_UCSC

for species in hg19 rheMac2 mm9 rheMac10 mm39 rn7 canFam6 equCab3
do
    wget https://hgdownload.cse.ucsc.edu/goldenpath/${species}/database/cpgIslandExt.txt.gz
    gunzip cpgIslandExt.txt.gz
    mv cpgIslandExt.txt ${species}_CGI_standardUCSC.bed
done

##### STEP 1:

# MAKE BED FILES FOR ALL SETS with names in 4th column indicating list of origin

mkdir /gpfs/gibbs/pi/noonan/ak2267/Revisions/Model-based_CGIs/20240305_overlaps_with_more_info
cd /gpfs/gibbs/pi/noonan/ak2267/Revisions/Model-based_CGIs/20240305_overlaps_with_more_info

mkdir Bed

# STANDARD UCSC FILES:
# in /home/ak2267/genomes/CGI/Standard_UCSC

for species in hg19 mm9 rheMac10 mm39 rn7 canFam6 equCab3
do
    cat /home/ak2267/genomes/CGI/Standard_UCSC/${species}_CGI_standardUCSC.bed | awk '{ print $2"\t"$3"\t"$4"\tUCSC_Standard" }' - > Bed/${species}_CGI_UCSC_Standard.bed
done

for species in rheMac2
do
    cat /home/ak2267/genomes/CGI/Standard_UCSC/${species}_CGI_standardUCSC.bed | awk '{ print $1"\t"$2"\t"$3"\tUCSC_Standard" }' - > Bed/${species}_CGI_UCSC_Standard.bed
done

# UCSC_AL FILES:
# in /home/ak2267/genomes/CGI/UCSC_AL/

for species in hg19 rheMac2 mm9 rheMac10 mm39 rn7 canFam6 equCab3
do
    cat /home/ak2267/genomes/CGI/UCSC_AL/${species}_CGIsAL_preMerge.bed | awk '{ print $1"\t"$2"\t"$3"\tUCSC_AL" }' - > Bed/${species}_CGI_UCSC_AL.bed
done

# MODEL-BASED FILES
# merged versions are in /gpfs/gibbs/pi/noonan/ak2267/Revisions/Model-based_CGIs/Simple_overlaps/Bed
# I have to re-make these to work with unmerged versions
# HERE I AM CONVERTING FILES FROM 1-BASED TO 0-BASED

mkdir Model-based_preLift
for species in hg19 rheMac2 mm9 mm10 rn4 canFam2 equCab2
do
    cut -f 1,2,3 /home/ak2267/genomes/CGI/Model_based/model-based-cpg-islands-${species}.txt | sed '1d' | \
    awk '{ print $1"\t"($2 - 1)"\t"$3"\tModel_based" }' > Model-based_preLift/${species}_model-based_CGIs.bed
done

# lift as necessary (see text under STEP 1)
module load liftOver/2023-05-23
# needed to get interactive node with more memory: salloc -t 6:00:00 --mem=16G

liftOver -minMatch=0.8 Model-based_preLift/rheMac2_model-based_CGIs.bed /home/ak2267/genomes/chain/rheMac2ToRheMac8.over.chain.gz Model-based_preLift/rheMac2_model-based_CGIs_LOrheMac8.bed unMapped # 
liftOver -minMatch=0.8 Model-based_preLift/rheMac2_model-based_CGIs_LOrheMac8.bed /home/ak2267/genomes/chain/rheMac8ToRheMac10.over.chain.gz Bed/rheMac10_CGI_Model_based.bed unMapped # 
liftOver -minMatch=0.8 Model-based_preLift/mm10_model-based_CGIs.bed /home/ak2267/genomes/chain/mm10ToMm39.over.chain.gz Bed/mm39_CGI_Model_based.bed unMapped # 
liftOver -minMatch=0.8 Model-based_preLift/rn4_model-based_CGIs.bed /home/ak2267/genomes/chain/rn4ToRn6.over.chain.gz Model-based_preLift/rn4_model-based_CGIs_LOrn6.bed unMapped # 
liftOver -minMatch=0.8 Model-based_preLift/rn4_model-based_CGIs_LOrn6.bed /home/ak2267/genomes/chain/rn6ToRn7.over.chain.gz Bed/rn7_CGI_Model_based.bed unMapped # 
liftOver -minMatch=0.8 Model-based_preLift/canFam2_model-based_CGIs.bed /home/ak2267/genomes/chain/canFam2ToCanFam6.over.chain.gz Bed/canFam6_CGI_Model_based.bed unMapped # 
liftOver -minMatch=0.8 Model-based_preLift/equCab2_model-based_CGIs.bed /home/ak2267/genomes/chain/equCab2ToEquCab3.over.chain.gz Bed/equCab3_CGI_Model_based.bed unMapped # 

# plus just move hg19, mm9, and rheMac2 over
for species in hg19 rheMac2 mm9
do
    cp Model-based_preLift/${species}_model-based_CGIs.bed Bed/${species}_CGI_Model_based.bed
done

# THIS MAKES FILES FOR EACH SPECIES IN /gpfs/gibbs/pi/noonan/ak2267/Revisions/Model-based_CGIs/20240305_overlaps_with_more_info/Bed/
# {species}_CGI_UCSC_Standard.bed     - has "UCSC_Standard" in col 4
# {species}_CGI_UCSC_AL.bed           - has "UCSC_AL" in col 4
# {species}_CGI_Model_based.bed       - has "Model_based" in col 4

##### STEP 2:

# RESTRICT TO CHROMOSOMES USED IN ALL ANALYSES:

for species in hg19 rheMac2 mm9 rheMac10 mm39 rn7 canFam6 equCab3
do 
    echo $species
    cut -f 1 Bed/${species}_CGI_UCSC_Standard.bed | sort -u | wc -l
    cut -f 1 Bed/${species}_CGI_UCSC_AL.bed | sort -u | wc -l
    cut -f 1 Bed/${species}_CGI_Model_based.bed | sort -u | wc -l
done

for species in mm39
do 
    echo $species
    echo '#### UCSC_Standard'
    cut -f 1 Bed/${species}_CGI_UCSC_Standard.bed | sort -u 
    echo '#### UCSC_AL'
    cut -f 1 Bed/${species}_CGI_UCSC_AL.bed | sort -u 
    echo '#### Model_based'
    cut -f 1 Bed/${species}_CGI_Model_based.bed | sort -u 
done

# Use python script to make bed files that only have CGIs on chromosomes present in all 3 files

mkdir Bed_commonChrom

for species in hg19 rheMac2 mm9 rheMac10 mm39 rn7 canFam6 equCab3
do 
    python ak20240305_restrict_to_common_chromosomes.py Bed/${species}_CGI_UCSC_Standard.bed Bed/${species}_CGI_UCSC_AL.bed Bed/${species}_CGI_Model_based.bed Bed_commonChrom/${species}_CGI_UCSC_Standard.bed Bed_commonChrom/${species}_CGI_UCSC_AL.bed Bed_commonChrom/${species}_CGI_Model_based.bed
done

for species in hg19 rheMac2 mm9 rheMac10 mm39 rn7 canFam6 equCab3
do 
    echo $species
    cut -f 1 Bed_commonChrom/${species}_CGI_UCSC_Standard.bed | sort -u | wc -l
    cut -f 1 Bed_commonChrom/${species}_CGI_UCSC_AL.bed | sort -u | wc -l
    cut -f 1 Bed_commonChrom/${species}_CGI_Model_based.bed | sort -u | wc -l
done

# MAKES FILES IN /gpfs/gibbs/pi/noonan/ak2267/Revisions/Model-based_CGIs/20240305_overlaps_with_more_info/Bed_commonChrom/
# THAT ONLY HAVE CGIs ON CHROMOSOMES PRESENT IN ALL 3 CGI LISTS
# {species}_CGI_UCSC_Standard.bed     
# {species}_CGI_UCSC_AL.bed           
# {species}_CGI_Model_based.bed       

##### STEP 3: 

# MERGE ACROSS ALL LISTS

for species in hg19 rheMac2 mm9 rheMac10 mm39 rn7 canFam6 equCab3
do 
    cat Bed_commonChrom/${species}_CGI_UCSC_Standard.bed Bed_commonChrom/${species}_CGI_UCSC_AL.bed Bed_commonChrom/${species}_CGI_Model_based.bed > ${species}_cat.bed
    sort -k1,1 -k2,2n ${species}_cat.bed > ${species}_sort.bed
    bedtools merge -c 4 -o distinct -i ${species}_sort.bed > ${species}_mergedCGIs.bed
done

rm *_cat.bed
rm *_sort.bed

species=hg19

grep Model_based ${species}_mergedCGIs.bed | grep -v UCSC_AL - | grep -v UCSC_Standard | wc -l
grep UCSC_AL ${species}_mergedCGIs.bed | grep -v Model_based - | grep -v UCSC_Standard | wc -l
grep UCSC_Standard ${species}_mergedCGIs.bed | grep -v Model_based - | grep -v UCSC_AL | wc -l

grep Model_based,UCSC_AL ${species}_mergedCGIs.bed | grep -v UCSC_Standard - | wc -l
grep Model_based,UCSC_Standard ${species}_mergedCGIs.bed | grep -v UCSC_AL - | wc -l
grep UCSC_AL,UCSC_Standard ${species}_mergedCGIs.bed | grep -v Model_based - | wc -l

grep Model_based,UCSC_AL,UCSC_Standard ${species}_mergedCGIs.bed | wc -l


##### STEP 4: 
# USE FACOUNT TO GET INFO ON INTERVALS: length, CpG num, GC_pct, CpG Obs/Exp
# DO NEGATIVE OVERLAP WITH FEATURES TO GET oCGIs
# Use python script to get results into tables for use in R

mkdir faCount

# get masked fasta files that I don't have
cd /gpfs/gibbs/pi/noonan/ak2267/genomes/
for genome in rheMac10 mm39 rn7 canFam6 equCab3
do
    echo 'cd /gpfs/gibbs/pi/noonan/ak2267/genomes/ ; wget https://hgdownload.soe.ucsc.edu/goldenPath/'${genome}'/bigZips/'${genome}'.fa.masked.gz ; gunzip '${genome}'.fa.masked.gz'
done >> 240305_downloadMaskedGenomes_jobFile.txt
dsq --job-file 240305_downloadMaskedGenomes_jobFile.txt --mem-per-cpu 5G -c 1 --mail-type FAIL,END
sbatch dsq-240305_downloadMaskedGenomes_jobFile-2024-03-05.sh # 22875514

# run bedtools getfasta and faCount - very quick on interactive node (<20 sec)
cd /gpfs/gibbs/pi/noonan/ak2267/Revisions/Model-based_CGIs/20240305_overlaps_with_more_info

for species in hg19 rheMac2 mm9 rheMac10 mm39 rn7 canFam6 equCab3
do
    echo $species
    bedtools getfasta -name -fi /gpfs/gibbs/pi/noonan/ak2267/genomes/${species}.fa.masked -bed ${species}_mergedCGIs.bed > faCount/${species}_mergedCGIs.fa ; faCount faCount/${species}_mergedCGIs.fa > faCount/${species}_mergedCGIs.faCount
done

# intersect merged files with feature annotations to identify oCGIs
mkdir noFeatures

for species in hg19 mm9 mm39
do
    echo $species
    bedtools intersect -v -a ${species}_mergedCGIs.bed \
	-b /home/ak2267/genomes/RefSeq/featureAnnotations/${species}_allFeatures.bed \
	-b /home/ak2267/genomes/FANTOM/FANTOM_TSS_${species}.bed \
	-b /home/ak2267/genomes/blacklist/blacklist_${species}.bed > noFeatures/${species}_mergedCGIs_noFeatures.bed
done

for species in rheMac2 rheMac10 rn7 canFam6 equCab3
do
    echo $species
    bedtools intersect -v -a ${species}_mergedCGIs.bed \
	 -b /home/ak2267/genomes/RefSeq/featureAnnotations/${species}_allFeatures.bed > noFeatures/${species}_mergedCGIs_noFeatures.bed
done


# MODIFICATIONS 19/3/24
# re-run getfasta on unmasked genomes (takes <20 sec on interactive node)
for species in hg19 rheMac2 mm9 rheMac10 mm39 rn7 canFam6 equCab3
do
    echo $species
    bedtools getfasta -name -fi /gpfs/gibbs/pi/noonan/ak2267/genomes/${species}.fa -bed ${species}_mergedCGIs.bed > faCount/${species}_mergedCGIs.fa ; faCount faCount/${species}_mergedCGIs.fa > faCount/${species}_mergedCGIs.faCount
done

# then modify python script to incorporate N counts
# use python script to make output tables for R - runs in < 30 sec per species
# a FEW lines have a 0 value for a C or a G - I think these are from liftover related errors (i.e. masking changes between genome versions). the numbers are:
# only happens ONCE in rn7 and ONCE in equCab3
grep NA CGI_comparisons_20240319/*
# CGI_comparisons_20240319/equCab3_CGI_comparisons.txt:chr28:18726714-18726804	1	0	0	90	90	0.0	0	NA	1
# CGI_comparisons_20240319/rn7_CGI_comparisons.txt:chr7:2455916-2456103	1	0	0	187	187	0.0	0	NA	1


# ALSO: incorporate intersection with rmsk (RepeatMasker) track - get number of bases overlapping a repeat for each CGI
mkdir Intersect_rmsk
for species in hg19 rheMac2 mm9 rheMac10 mm39 rn7 canFam6 equCab3
do
    echo $species
    bedtools intersect -wao -a ${species}_mergedCGIs.bed -b /home/ak2267/genomes/rmsk/rmsk_hg19_merged.bed > Intersect_rmsk/${species}_mergedCGIs_intersectRmsk.txt
done


# runs in ~30 sec per species - gets longer with rmsk step
for species in hg19 rheMac2 mm9 rheMac10 mm39 rn7 canFam6 equCab3
do
    echo $species
    python ak20240319_makeTablesForR.py faCount/${species}_mergedCGIs.faCount noFeatures/${species}_mergedCGIs_noFeatures.bed Intersect_rmsk/${species}_mergedCGIs_intersectRmsk.txt CGI_comparisons_20240319/${species}_CGI_comparisons.txt
done

# download for use in R
tar -zcvf CGI_comparisons_20240319.gz CGI_comparisons_20240319/

cd Desktop/Yale/\!Writing/Kocher_CpG/GB_Revisions/Revised_code/Step12_revisions/Model-based_CGIs 
scp ak2267@mccleary.ycrc.yale.edu:/gpfs/gibbs/pi/noonan/ak2267/Revisions/Model-based_CGIs/20240305_overlaps_with_more_info/CGI_comparisons_20240319.gz .

