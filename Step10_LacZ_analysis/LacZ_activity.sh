# 7/15/22
# Purpose: take enhancers with each mark (H3K4me3, H3K27ac, H3K4me1, H3K4me2) in embryonic mouse tissues (e11.5)
# See how many are also active in LacZ experiments, as a way to demonstrate that the marks predict enhancer activity

####### STEP 1: acquire entire VISTA database from Len Pennacchio (see 7/15/22 email) - THIS IS IN mm10!
# VISTA:
# https://enhancer.lbl.gov/
# ENCODE paper that tested 151 enhancers at e11.5, selected based on DNase and H3K27ac signal (3 tiers) in midbrain, hindbrain, and limb:
# https://www.nature.com/articles/s41586-020-2493-4#Sec12
# ENCODE paper that tested 150 enhancers at e12.5, selected based on H3K27ac signal (3 tiers) in forebrain, heart, and limb:
# https://www.nature.com/articles/s41586-020-2093-3#Sec6
# It looks like both ENCODE papers deposited their tested enhancers into VISTA with names accordingly (mm####)

# copy files to Ruddle
cd /Users/acadiak/Desktop/CGI/ENCODE_LacZ
scp acadia/2022-07-15_VISTA.tsv ak2267@ruddle.hpc.yale.edu:/gpfs/gibbs/pi/noonan/ak2267/VISTA
scp acadia/tissue.dictionary.csv ak2267@ruddle.hpc.yale.edu:/gpfs/gibbs/pi/noonan/ak2267/VISTA

# to Ruddle
cd /gpfs/gibbs/pi/noonan/ak2267/VISTA

# script to process into bed files for each species x timePoint x backbone, then delete empty files since I didn't do so from within the python script
python makeManyBedsFromVISTA.py 2022-07-15_VISTA.tsv tissue.dictionary.csv
find . -empty -delete 

####### STEP 2: download bed files for each mark in e11.5 tissues - THIS IS IN mm10!
# {e11.5, e12.5} x {forebrain, midbrain, hindbrain, heart, limb, embryonic facial prominence, neural tube, liver} x {H3K4me3, H3K27ac, H3K4me1, H3K4me2}

cd /gpfs/gibbs/pi/noonan/ak2267/VISTA/ENCODE

# start with e11.5, restrict to {forebrain, midbrain, hindbrain, heart, limb} {H3K4me3, H3K27ac, H3K4me1, H3K4me2}
# H3K4me3
wget https://www.encodeproject.org/files/ENCFF958GPD/@@download/ENCFF958GPD.bed.gz ; mv ENCFF958GPD.bed.gz e11.5_forebrain_H3K4me3.bed.gz
wget https://www.encodeproject.org/files/ENCFF651IYR/@@download/ENCFF651IYR.bed.gz ; mv ENCFF651IYR.bed.gz e11.5_midbrain_H3K4me3.bed.gz
wget https://www.encodeproject.org/files/ENCFF445DQR/@@download/ENCFF445DQR.bed.gz ; mv ENCFF445DQR.bed.gz e11.5_hindbrain_H3K4me3.bed.gz
wget https://www.encodeproject.org/files/ENCFF908XCE/@@download/ENCFF908XCE.bed.gz ; mv ENCFF908XCE.bed.gz e11.5_heart_H3K4me3.bed.gz
wget https://www.encodeproject.org/files/ENCFF928TGM/@@download/ENCFF928TGM.bed.gz ; mv ENCFF928TGM.bed.gz e11.5_limb_H3K4me3.bed.gz

# H3K27ac
wget https://www.encodeproject.org/files/ENCFF897EEM/@@download/ENCFF897EEM.bed.gz ; mv ENCFF897EEM.bed.gz e11.5_forebrain_H3K27ac.bed.gz
wget https://www.encodeproject.org/files/ENCFF566DFK/@@download/ENCFF566DFK.bed.gz ; mv ENCFF566DFK.bed.gz e11.5_midbrain_H3K27ac.bed.gz
wget https://www.encodeproject.org/files/ENCFF203QTV/@@download/ENCFF203QTV.bed.gz ; mv ENCFF203QTV.bed.gz e11.5_hindbrain_H3K27ac.bed.gz
wget https://www.encodeproject.org/files/ENCFF236UMU/@@download/ENCFF236UMU.bed.gz ; mv ENCFF236UMU.bed.gz e11.5_heart_H3K27ac.bed.gz
wget https://www.encodeproject.org/files/ENCFF016BEF/@@download/ENCFF016BEF.bed.gz ; mv ENCFF016BEF.bed.gz e11.5_limb_H3K27ac.bed.gz

# H3K4me1
wget https://www.encodeproject.org/files/ENCFF147OKD/@@download/ENCFF147OKD.bed.gz ; mv ENCFF147OKD.bed.gz e11.5_forebrain_H3K4me1.bed.gz
wget https://www.encodeproject.org/files/ENCFF202HIO/@@download/ENCFF202HIO.bed.gz ; mv ENCFF202HIO.bed.gz e11.5_midbrain_H3K4me1.bed.gz
wget https://www.encodeproject.org/files/ENCFF098IGX/@@download/ENCFF098IGX.bed.gz ; mv ENCFF098IGX.bed.gz e11.5_hindbrain_H3K4me1.bed.gz
wget https://www.encodeproject.org/files/ENCFF218FKJ/@@download/ENCFF218FKJ.bed.gz ; mv ENCFF218FKJ.bed.gz e11.5_heart_H3K4me1.bed.gz
wget https://www.encodeproject.org/files/ENCFF255WOG/@@download/ENCFF255WOG.bed.gz ; mv ENCFF255WOG.bed.gz e11.5_limb_H3K4me1.bed.gz

# H3K4me2
wget https://www.encodeproject.org/files/ENCFF047OVD/@@download/ENCFF047OVD.bed.gz ; mv ENCFF047OVD.bed.gz e11.5_forebrain_H3K4me2.bed.gz
wget https://www.encodeproject.org/files/ENCFF132QFU/@@download/ENCFF132QFU.bed.gz ; mv ENCFF132QFU.bed.gz e11.5_midbrain_H3K4me2.bed.gz
wget https://www.encodeproject.org/files/ENCFF835VQG/@@download/ENCFF835VQG.bed.gz ; mv ENCFF835VQG.bed.gz e11.5_hindbrain_H3K4me2.bed.gz
wget https://www.encodeproject.org/files/ENCFF316YSJ/@@download/ENCFF316YSJ.bed.gz ; mv ENCFF316YSJ.bed.gz e11.5_heart_H3K4me2.bed.gz
wget https://www.encodeproject.org/files/ENCFF703AUU/@@download/ENCFF703AUU.bed.gz ; mv ENCFF703AUU.bed.gz e11.5_limb_H3K4me2.bed.gz


####### STEP 3:  intersect VISTA with intronic/intergenic peaks and CGIs

# LacZ file I'm using here (generated in step 1 above):
# mm10_e11.5_Hsp68_LacZ.bed

# CGI file was re-made 4/21/23 using text involving temp_mm39_LOmm10.bed in 220717_UCSC_AL_CpGislands.sh
# in order to have MERGED CGI tracks (recent update to pipeline)

cd /gpfs/gibbs/pi/noonan/ak2267/VISTA

mkdir intersect_VISTA-centric

for tissue in forebrain midbrain hindbrain heart limb
do
    for mark in H3K4me3 H3K27ac H3K4me1 H3K4me2
    do
	echo $tissue $mark
	bedtools intersect -v -a mm10_e11.5_Hsp68_LacZ.bed -b AK_only/mm39_allFeatures_LOmm10.bed | bedtools intersect -wao -a - -b ENCODE/e11.5_${tissue}_${mark}.bed > intersect_VISTA-centric/e11.5_${tissue}_${mark}_intersectPeaks.txt
	bedtools intersect -v -a mm10_e11.5_Hsp68_LacZ.bed -b AK_only/mm39_allFeatures_LOmm10.bed | bedtools intersect -wo -a - -b /home/ak2267/genomes/CGI/UCSC_AL/temp_mm39_LOmm10.bed > intersect_VISTA-centric/e11.5_${tissue}_${mark}_intersectCGIs.txt
    done;done

######## STEP 4: python script to summarize results, but VISTA-centric instead of peak-centric
# The reason for this is to capture a "background rate" using VISTA elements that do NOT overlap peaks in the tissue of interest

# goal: table for each tissue x mark = 5 x 4 = 20
# then can manipulate in R to make barplots with different overlap cutoffs
# and can do Fisher's exact tests for enrichment compared to elements not overlapping peaks

mkdir summaryTables_VISTA-centric

for tissue in forebrain midbrain hindbrain heart limb
do
    for mark in H3K4me3 H3K27ac H3K4me1 H3K4me2
    do
	echo $tissue $mark
	python makeSummaryTableForR_VISTA-centric.py ${tissue} ${mark} > summaryTables_VISTA-centric/e11.5_${tissue}_${mark}_summaryTable.txt
    done;done
    
####### STEP 5

tar -zcvf summaryTables_VISTA-centric.gz summaryTables_VISTA-centric

# to laptop
cd /Users/acadiak/Desktop/CGI/ENCODE_LacZ
scp ak2267@ruddle.hpc.yale.edu:/gpfs/gibbs/pi/noonan/ak2267/VISTA/summaryTables_VISTA-centric.gz .

# then analyze in R





