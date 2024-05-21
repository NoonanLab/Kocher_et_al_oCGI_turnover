# 2/27/2024
# Purpose: determine overlap of oCGIs with ENCODE cCRE categories
# and generate bar plots showing this
# This is revised manuscript Fig S6

# Make directory
mkdir /gpfs/gibbs/pi/noonan/ak2267/Revisions
mkdir /gpfs/gibbs/pi/noonan/ak2267/Revisions/ENCODE_cCREs

# Download ENCODE cCRE annotations in mm10 and hg38 and convert from bigBed to bed
cd /gpfs/gibbs/pi/noonan/ak2267/Revisions/ENCODE_cCREs
module load Kent_tools/411-GCC-10.2.0
bigBedToBed http://hgdownload.soe.ucsc.edu/gbdb/mm10/encode3/ccre/encodeCcreCombined.bb stdout > encodeCcreCombined_mm10.bed
bigBedToBed http://hgdownload.soe.ucsc.edu/gbdb/hg38/encode3/ccre/encodeCcreCombined.bb stdout > encodeCcreCombined_hg38.bed

# Make files with name of element type in 4th column
cut -f 1,2,3,10 encodeCcreCombined_hg38.bed > encodeCcreCombined_hg38_4col.bed
cut -f 1,2,3,10 encodeCcreCombined_mm10.bed > encodeCcreCombined_mm10_4col.bed

# Lift mouse from mm10 to mm39
liftOver -minMatch=0.8 encodeCcreCombined_mm10_4col.bed /home/ak2267/genomes/chain/mm10ToMm39.over.chain.gz \
encodeCcreCombined_mm10_LOmm39.bed unMapped
# 7834 / 343731 fail to map

# Generate oCGI lists by removing CGIs overlapping annotated features
for species in hg38 mm39
do
bedtools intersect -v -a /home/ak2267/genomes/CGI/UCSC_AL/${species}_CGIsAL.bed \
	-b /home/ak2267/genomes/RefSeq/featureAnnotations/${species}_allFeatures.bed \
	-b /home/ak2267/genomes/FANTOM/FANTOM_TSS_${species}.bed \
	-b /home/ak2267/genomes/blacklist/blacklist_${species}.bed \
	| awk '{ print $1"\t"$2"\t"$3"\tCGI_"NR }' - \
	> ${species}_CGIsAL_noFeatures.bed
done

# Overlap oCGIs with cCREs
bedtools intersect -wa -wb -a mm39_CGIsAL_noFeatures.bed -b encodeCcreCombined_mm10_LOmm39.bed | wc -l
# 7052
bedtools intersect -u -a mm39_CGIsAL_noFeatures.bed -b encodeCcreCombined_mm10_LOmm39.bed | wc -l
# 5732
# Here we can see that ~1300 oCGIs overlap multiple cCREs

bedtools intersect -wa -wb -a hg38_CGIsAL_noFeatures.bed -b encodeCcreCombined_hg38_4col.bed | wc -l
# 12190
bedtools intersect -u -a hg38_CGIsAL_noFeatures.bed -b encodeCcreCombined_hg38_4col.bed | wc -l
# 9360
# Here we can see that ~3000 oCGIs overlap multiple cCREs

# DECISION: report ALL cCREs that come out of these overlaps, and make stacked bars of their identities


###### Count numbers of each cCRE type

# MOUSE
for cCRE in CTCF-only dELS DNase-H3K4me3 pELS PLS
do
    shortName=$(bedtools intersect -wa -wb -a mm39_CGIsAL_noFeatures.bed -b encodeCcreCombined_mm10_LOmm39.bed | grep ${cCRE} | wc -l)
    longName=$(bedtools intersect -wa -wb -a mm39_CGIsAL_noFeatures.bed -b encodeCcreCombined_mm10_LOmm39.bed | grep ${cCRE},CTCF-bound | wc -l)
    difference=$(expr $shortName - $longName)
    
    echo ${cCRE} $difference
    echo ${cCRE},CTCF-bound $longName
done

# CTCF-only 0
# CTCF-only,CTCF-bound 340
# dELS 4385
# dELS,CTCF-bound 1247
# DNase-H3K4me3 254
# DNase-H3K4me3,CTCF-bound 229
# pELS 430
# pELS,CTCF-bound 121
# PLS 35
# PLS,CTCF-bound 11

# HUMAN
for cCRE in CTCF-only dELS DNase-H3K4me3 pELS PLS
do
    shortName=$(bedtools intersect -wa -wb -a hg38_CGIsAL_noFeatures.bed -b encodeCcreCombined_hg38_4col.bed | grep ${cCRE} | wc -l)
    longName=$(bedtools intersect -wa -wb -a hg38_CGIsAL_noFeatures.bed -b encodeCcreCombined_hg38_4col.bed | grep ${cCRE},CTCF-bound | wc -l)
    difference=$(expr $shortName - $longName)
    
    echo ${cCRE} $difference
    echo ${cCRE},CTCF-bound $longName
done

# CTCF-only 0
# CTCF-only,CTCF-bound 390
# dELS 4932
# dELS,CTCF-bound 4913
# DNase-H3K4me3 487
# DNase-H3K4me3,CTCF-bound 510
# pELS 398
# pELS,CTCF-bound 458
# PLS 46
# PLS,CTCF-bound 56

###### Count totals of each cCRE type in my intronic-intergenic genome space

# MOUSE

# get intronic-intergenic cCREs
bedtools intersect -v -a encodeCcreCombined_mm10_LOmm39.bed \
	-b /home/ak2267/genomes/RefSeq/featureAnnotations/mm39_allFeatures.bed \
	-b /home/ak2267/genomes/FANTOM/FANTOM_TSS_mm39.bed \
	-b /home/ak2267/genomes/blacklist/blacklist_mm39.bed \
	> mm39_cCREs_noFeatures.bed
	# 194003 / 339814 make it through

for cCRE in CTCF-only dELS DNase-H3K4me3 pELS PLS
do
    shortName=$(grep ${cCRE} mm39_cCREs_noFeatures.bed | wc -l)
    longName=$(grep ${cCRE},CTCF-bound mm39_cCREs_noFeatures.bed | wc -l)
    difference=$(expr $shortName - $longName)
    
    echo ${cCRE} $difference
    echo ${cCRE},CTCF-bound $longName
done

# CTCF-only 0
# CTCF-only,CTCF-bound 19332
# dELS 132191
# dELS,CTCF-bound 18671
# DNase-H3K4me3 3719
# DNase-H3K4me3,CTCF-bound 1840
# pELS 14873
# pELS,CTCF-bound 2916
# PLS 364
# PLS,CTCF-bound 97

# HUMAN

# get intronic-intergenic cCREs
bedtools intersect -v -a encodeCcreCombined_hg38_4col.bed \
	-b /home/ak2267/genomes/RefSeq/featureAnnotations/hg38_allFeatures.bed \
	-b /home/ak2267/genomes/FANTOM/FANTOM_TSS_hg38.bed \
	-b /home/ak2267/genomes/blacklist/blacklist_hg38.bed \
	> hg38_cCREs_noFeatures.bed
	# 617674 / 926535 make it through
	
for cCRE in CTCF-only dELS DNase-H3K4me3 pELS PLS
do
    shortName=$(grep ${cCRE} hg38_cCREs_noFeatures.bed | wc -l)
    longName=$(grep ${cCRE},CTCF-bound hg38_cCREs_noFeatures.bed | wc -l)
    difference=$(expr $shortName - $longName)
    
    echo ${cCRE} $difference
    echo ${cCRE},CTCF-bound $longName
done

# CTCF-only 0
# CTCF-only,CTCF-bound 49651
# dELS 352500
# dELS,CTCF-bound 155287
# DNase-H3K4me3 13586
# DNase-H3K4me3,CTCF-bound 6703
# pELS 22018
# pELS,CTCF-bound 16700
# PLS 620
# PLS,CTCF-bound 609

###### Count totals across entire cCRE datasets

# MOUSE

for cCRE in CTCF-only dELS DNase-H3K4me3 pELS PLS
do
    shortName=$(grep ${cCRE} encodeCcreCombined_mm10_LOmm39.bed | wc -l)
    longName=$(grep ${cCRE},CTCF-bound encodeCcreCombined_mm10_LOmm39.bed | wc -l)
    difference=$(expr $shortName - $longName)
    
    echo ${cCRE} $difference
    echo ${cCRE},CTCF-bound $longName
done

# CTCF-only 0
# CTCF-only,CTCF-bound 23836
# dELS 181172
# dELS,CTCF-bound 27867
# DNase-H3K4me3 6982
# DNase-H3K4me3,CTCF-bound 3401
# pELS 56411
# pELS,CTCF-bound 16383
# PLS 13855
# PLS,CTCF-bound 9907

# HUMAN

for cCRE in CTCF-only dELS DNase-H3K4me3 pELS PLS
do
    shortName=$(grep ${cCRE} encodeCcreCombined_hg38_4col.bed | wc -l)
    longName=$(grep ${cCRE},CTCF-bound encodeCcreCombined_hg38_4col.bed | wc -l)
    difference=$(expr $shortName - $longName)
    
    echo ${cCRE} $difference
    echo ${cCRE},CTCF-bound $longName
done

# CTCF-only 0
# CTCF-only,CTCF-bound 56766
# dELS 448981
# dELS,CTCF-bound 218618
# DNase-H3K4me3 16737
# DNase-H3K4me3,CTCF-bound 8800
# pELS 64421
# pELS,CTCF-bound 77409
# PLS 7582
# PLS,CTCF-bound 27221

# Copy results into Excel and then export as text file for use in R
# Excel file: oCGI_cCRE_overlap.xlsx
# Text file: oCGI_cCRE_overlap.txt
# R script with plots: oCGI_cCRE_overlap.R
