# 3/9/2024
# Purpose: take all mouse cCREs 
# Separate by: 1) overlapping a histone mod (H3K4me3, H3K27ac, H3K4me1, H3K4me2) in each of 5 e11.5 tissues I used before
#              2) overlapping an oCGI
# Count promoter interactions

# THIS WILL BE REVISED FIG S4

# work here:
cd /gpfs/gibbs/pi/noonan/ak2267/Revisions/Looping

# Count loops in each case
# This is an orthogonal approach to the LacZ reporter results I used before
# Many of these files were generated previously in LacZ_activity.sh
# For Fig S3

####### STEP 0: get mouse cCREs & restrict to those in my "oCGI space"
# by removing those overlapping features & masked regions

# cCREs are here:
/gpfs/gibbs/pi/noonan/ak2267/Revisions/ENCODE_cCREs/encodeCcreCombined_mm10.bed

# remove those overlapping features and blacklist space
bedtools intersect -v -a /gpfs/gibbs/pi/noonan/ak2267/Revisions/ENCODE_cCREs/encodeCcreCombined_mm10.bed \
    -b /gpfs/gibbs/pi/noonan/ak2267/VISTA/AK_only/mm39_allFeatures_LOmm10.bed \
	-b /home/ak2267/genomes/FANTOM/TSS_mouse_mm10.bed \
	-b /home/ak2267/genomes/blacklist/blacklist_mm10.bed | awk '{ print $1"\t"$2"\t"$3"\t"$4";"$10 }' |\
	grep -v fix | grep -v alt > encodeCcreCombined_mm10_noFeatures.bed
	
# use this file:
encodeCcreCombined_mm10_noFeatures.bed

####### STEP 1: LOOPS INSTEAD OF VISTA

# download from E-MTAB-2414: https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-2414
wget https://www.ebi.ac.uk/biostudies/files/E-MTAB-2414/FLC_promoter_other_significant_interactions.txt

# turn into bed file with coordinates of the *interaction* (not the promoter) and gene ENSG in column 4
# some have not_promoter in column_4
sed 's/\"//g' FLC_promoter_other_significant_interactions.txt | awk '{ print $7"\t"$8"\t"$9"\t"$5 }' | sed 1d > FLC_promoter_other_significant_interactions_mm9.bed
liftOver -minMatch=0.8 FLC_promoter_other_significant_interactions_mm9.bed /home/ak2267/genomes/chain/mm9ToMm10.over.chain.gz FLC_promoter_other_significant_interactions_mm10.bed unMapped

# use this file:
/gpfs/gibbs/pi/noonan/ak2267/Revisions/Looping/FLC_promoter_other_significant_interactions_mm10.bed

####### STEP 2: HISTONE MOD DATA in mm10
# {e14.5 to match HiC data} x {liver} x {H3K4me3, H3K27ac, H3K4me1, H3K4me2}

# download to here (previously I was working with e11.5 data to match the LacZ readout):
cd /gpfs/gibbs/pi/noonan/ak2267/VISTA/ENCODE
wget https://www.encodeproject.org/files/ENCFF384FJW/@@download/ENCFF384FJW.bed.gz ; mv ENCFF384FJW.bed.gz e14.5_liver_H3K27ac.bed.gz
wget https://www.encodeproject.org/files/ENCFF510NON/@@download/ENCFF510NON.bed.gz ; mv ENCFF510NON.bed.gz e14.5_liver_H3K4me3.bed.gz
wget https://www.encodeproject.org/files/ENCFF627CYS/@@download/ENCFF627CYS.bed.gz ; mv ENCFF627CYS.bed.gz e14.5_liver_H3K4me2.bed.gz
wget https://www.encodeproject.org/files/ENCFF880DXF/@@download/ENCFF880DXF.bed.gz ; mv ENCFF880DXF.bed.gz e14.5_liver_H3K4me1.bed.gz
gunzip *.gz

####### STEP 3:  intersect cCREs with peaks, CGIs, and interactions

# CGI file:
/home/ak2267/genomes/CGI/UCSC_AL/temp_mm39_LOmm10.bed

cd /gpfs/gibbs/pi/noonan/ak2267/Revisions/Looping

mkdir intersect_peaks
mkdir intersect_CGIs
mkdir intersect_interactions

# intersect with peaks
for mark in H3K4me3 H3K27ac H3K4me1 H3K4me2
do
	echo $mark
	cut -f 1,2,3,4 /gpfs/gibbs/pi/noonan/ak2267/VISTA/ENCODE/e14.5_liver_${mark}.bed | bedtools intersect -wao -a encodeCcreCombined_mm10_noFeatures.bed -b - > intersect_peaks/e14.5_liver_${mark}_intersectPeaks.txt
done

# intersect with CGIs
bedtools intersect -wao -a encodeCcreCombined_mm10_noFeatures.bed -b /home/ak2267/genomes/CGI/UCSC_AL/temp_mm39_LOmm10.bed > intersect_CGIs/encodeCcreCombined_mm10_noFeatures_intersectCGIs.txt

# intersect with interactions
bedtools intersect -wao -a encodeCcreCombined_mm10_noFeatures.bed -b FLC_promoter_other_significant_interactions_mm10.bed > intersect_interactions/encodeCcreCombined_mm10_noFeatures_intersectInteractions.txt


######## STEP 4: python script to summarize results

# goal: table for each mark = 4
# then can manipulate in R to make barplots with different overlap cutoffs
# and can do Fisher's exact tests for enrichment compared to elements not overlapping peaks

mkdir summaryTables

for mark in H3K4me3 H3K27ac H3K4me1 H3K4me2
do
	echo $mark
	python makeSummaryTableForR_looping.py ${mark} > summaryTables/${mark}_summaryTable_looping.txt
done


####### STEP 5

tar -zcvf summaryTables_looping.gz summaryTables

# to laptop
#cd /Users/acadiak/Desktop/Yale/!Writing/Kocher_CpG/GB_Revisions/Revised_code/Step12_revisions/Loops
#scp ak2267@mccleary.ycrc.yale.edu:/gpfs/gibbs/pi/noonan/ak2267/Revisions/Looping/summaryTables_looping.gz .

# then analyze in R

# USE THIS VERSION FOR FIGURE #
# make cCREs larger ONLY FOR INTERACTION STEP and redo intersection with interactions
bedtools slop -b 1000 -i encodeCcreCombined_mm10_noFeatures.bed \
-g /home/ak2267/genomes/chrom.sizes/mm10.chrom.sizes | bedtools intersect -wao -a - \
-b FLC_promoter_other_significant_interactions_mm10.bed > intersect_interactions/encodeCcreCombined_mm10_noFeatures_intersectInteractions.txt

# new summary tables
mkdir summaryTables_pad1kb
for mark in H3K4me3 H3K27ac H3K4me1 H3K4me2
do
	echo $mark
	python makeSummaryTableForR_looping.py ${mark} > summaryTables_pad1kb/${mark}_summaryTable_looping.txt
done

tar -zcvf summaryTables_pad1kb.gz summaryTables_pad1kb

# to laptop
cd /Users/acadiak/Desktop/Yale/!Writing/Kocher_CpG/GB_Revisions/Revised_code/Step12_revisions/Loops
scp ak2267@mccleary.ycrc.yale.edu:/gpfs/gibbs/pi/noonan/ak2267/Revisions/Looping/summaryTables_pad1kb.gz .

# then analyze in R
