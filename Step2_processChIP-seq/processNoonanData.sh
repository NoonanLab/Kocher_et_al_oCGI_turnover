# 7/28/22
# Purpose: download and process Noonan lab ChIP data (developing brain from Reilly, Yin et al 2015 and developing limb from Cotney, Leng et al 2013)

############################## BRAIN (CORTEX) ######################################

# download peak files from Reilly, Yin 2015 (in hg19, rheMac2, mm9)
# H3K27ac
cd /home/ak2267/project/EnhancerClasses/hg19/ac
for timePoint in CS16 CS23 F2F F2O; do wget http://noonan.ycga.yale.edu/noonan_public/reilly2015/Track_Hub_ac/hg19/merge_${timePoint}_overlap.bb; done
cd /home/ak2267/project/EnhancerClasses/rheMac2/ac
for timePoint in e79F e79O; do wget http://noonan.ycga.yale.edu/noonan_public/reilly2015/Track_Hub_ac/rheMac2/merge_${timePoint}_overlap_rh.bb; done
cd /home/ak2267/project/EnhancerClasses/mm9/ac
for timePoint in 17F 17O e11 e14; do wget http://noonan.ycga.yale.edu/noonan_public/reilly2015/Track_Hub_ac/mm9/merge_${timePoint}_overlap_mm.bb; done
# H3K4me2
cd /home/ak2267/project/EnhancerClasses/hg19/me2
for timePoint in CS16 CS23 F2F F2O; do wget http://noonan.ycga.yale.edu/noonan_public/reilly2015/Track_Hub_me2/hg19/merge_${timePoint}_me2_overlap.bb; done
cd /home/ak2267/project/EnhancerClasses/rheMac2/me2
for timePoint in e79F e79O; do wget http://noonan.ycga.yale.edu/noonan_public/reilly2015/Track_Hub_me2/rheMac2/merge_${timePoint}_overlap_me2_rh.bb; done
cd /home/ak2267/project/EnhancerClasses/mm9/me2
for timePoint in 17F 17O e11 e14; do wget http://noonan.ycga.yale.edu/noonan_public/reilly2015/Track_Hub_me2/mm9/merge_${timePoint}_overlap_me2_mm.bb; done

# convert bigBed to bed
for mark in ac me2
do
    for file in $(ls hg19/${mark}/*.bb); do name=`basename $file .bb`; bigBedToBed ${file} hg19/${mark}/${name}.bed; done
    for file in $(ls rheMac2/${mark}/*.bb); do name=`basename $file .bb`; bigBedToBed ${file} rheMac2/${mark}/${name}.bed; done
    for file in $(ls mm9/${mark}/*.bb); do name=`basename $file .bb`; bigBedToBed ${file} mm9/${mark}/${name}.bed; done; done				 
# remove bigBed files
for directory in hg19 rheMac2 mm9; do rm ${directory}/*/*.bb; done

# get rhesus e55 time point (single replicate only), convert bigBed to bed, and rename in line with other files
cd /home/ak2267/project/EnhancerClasses/rheMac2/ac
wget http://noonan.ycga.yale.edu/noonan_public/reilly2015/Track_Hub_ac/rheMac2/e55_ac_rep1_regions.bb; bigBedToBed e55_ac_rep1_regions.bb merge_e55_overlap_rh.bed
cd /home/ak2267/project/EnhancerClasses/rheMac2/me2
wget http://noonan.ycga.yale.edu/noonan_public/reilly2015/Track_Hub_me2/rheMac2/e55_me2_rep1_regions.bb; bigBedToBed e55_me2_rep1_regions.bb merge_e55_overlap_me2_rh.bed

# add peak number to 4th column of peak files to keep track of which ones lift
for species in hg19 rheMac2 mm9
do
    for mark in ac me2
    do
	for file in ${species}/${mark}/*.bed
	do name=`basename $file .bed`; awk '{ print $1"\t"$2"\t"$3"\tPeak_"NR}' ${file} > ${species}/${mark}/${name}_named.bed
	done;done;done


# DOWNLOAD bigWigs for hg19, rheMac2, and mm9 - will need to quantify later using consensus peaks and bigWigAverageOverbed
mkdir /gpfs/gibbs/pi/noonan/ak2267/NoonanBrain
cd /gpfs/gibbs/pi/noonan/ak2267/NoonanBrain
mkdir hg19; mkdir rheMac2; mkdir mm9

# hg19
cd /gpfs/gibbs/pi/noonan/ak2267/NoonanBrain/hg19; mkdir ac; mkdir me2
for timePoint in CS16 CS23 F2F F2O
do
    cd /gpfs/gibbs/pi/noonan/ak2267/NoonanBrain/hg19/ac
    wget http://noonan.ycga.yale.edu/noonan_public/reilly2015/Track_Hub_ac/hg19/${timePoint}_k27ac_rep1.bw; mv ${timePoint}_k27ac_rep1.bw ${timePoint}_ac_rep1.bw
    wget http://noonan.ycga.yale.edu/noonan_public/reilly2015/Track_Hub_ac/hg19/${timePoint}_k27ac_rep2.bw; mv ${timePoint}_k27ac_rep2.bw ${timePoint}_ac_rep2.bw
    cd /gpfs/gibbs/pi/noonan/ak2267/NoonanBrain/hg19/me2
    wget http://noonan.ycga.yale.edu/noonan_public/reilly2015/Track_Hub_me2/hg19/${timePoint}_k4me2_rep1.bw; mv ${timePoint}_k4me2_rep1.bw ${timePoint}_me2_rep1.bw
    wget http://noonan.ycga.yale.edu/noonan_public/reilly2015/Track_Hub_me2/hg19/${timePoint}_k4me2_rep2.bw; mv ${timePoint}_k4me2_rep2.bw ${timePoint}_me2_rep2.bw
done

# rheMac2
cd /gpfs/gibbs/pi/noonan/ak2267/NoonanBrain/rheMac2; mkdir ac; mkdir me2
cd /gpfs/gibbs/pi/noonan/ak2267/NoonanBrain/rheMac2/ac
wget http://noonan.ycga.yale.edu/noonan_public/reilly2015/Track_Hub_ac/rheMac2/e55_k27ac_rep1.bw; mv e55_k27ac_rep1.bw e55_ac_rep1.bw
wget http://noonan.ycga.yale.edu/noonan_public/reilly2015/Track_Hub_ac/rheMac2/e79F_k27ac_rep1.bw; mv e79F_k27ac_rep1.bw e79F_ac_rep1.bw
wget http://noonan.ycga.yale.edu/noonan_public/reilly2015/Track_Hub_ac/rheMac2/e79F_k27ac_rep2.bw; mv e79F_k27ac_rep2.bw e79F_ac_rep2.bw
wget http://noonan.ycga.yale.edu/noonan_public/reilly2015/Track_Hub_ac/rheMac2/e79O_k27ac_rep1.bw; mv e79O_k27ac_rep1.bw e79O_ac_rep1.bw
wget http://noonan.ycga.yale.edu/noonan_public/reilly2015/Track_Hub_ac/rheMac2/e79O_k27ac_rep2.bw; mv e79O_k27ac_rep2.bw e79O_ac_rep2.bw
cd /gpfs/gibbs/pi/noonan/ak2267/NoonanBrain/rheMac2/me2
wget http://noonan.ycga.yale.edu/noonan_public/reilly2015/Track_Hub_me2/rheMac2/e55_k4me2_rep1.bw; mv e55_k4me2_rep1.bw e55_me2_rep1.bw
wget http://noonan.ycga.yale.edu/noonan_public/reilly2015/Track_Hub_me2/rheMac2/e79F_k4me2_rep1.bw; mv e79F_k4me2_rep1.bw e79F_me2_rep1.bw
wget http://noonan.ycga.yale.edu/noonan_public/reilly2015/Track_Hub_me2/rheMac2/e79F_k4me2_rep2.bw; mv e79F_k4me2_rep2.bw e79F_me2_rep2.bw
wget http://noonan.ycga.yale.edu/noonan_public/reilly2015/Track_Hub_me2/rheMac2/e79O_k4me2_rep1.bw; mv e79O_k4me2_rep1.bw e79O_me2_rep1.bw
wget http://noonan.ycga.yale.edu/noonan_public/reilly2015/Track_Hub_me2/rheMac2/e79O_k4me2_rep2.bw; mv e79O_k4me2_rep2.bw e79O_me2_rep2.bw

# mm9
cd /gpfs/gibbs/pi/noonan/ak2267/NoonanBrain/mm9; mkdir ac; mkdir me2
for timePoint in e11 e14 e17F e17O
do
    cd /gpfs/gibbs/pi/noonan/ak2267/NoonanBrain/mm9/ac
    wget http://noonan.ycga.yale.edu/noonan_public/reilly2015/Track_Hub_ac/mm9/${timePoint}_k27ac_rep1.bw; mv ${timePoint}_k27ac_rep1.bw ${timePoint}_ac_rep1.bw
    wget http://noonan.ycga.yale.edu/noonan_public/reilly2015/Track_Hub_ac/mm9/${timePoint}_k27ac_rep2.bw; mv ${timePoint}_k27ac_rep2.bw ${timePoint}_ac_rep2.bw
    cd /gpfs/gibbs/pi/noonan/ak2267/NoonanBrain/mm9/me2
    wget http://noonan.ycga.yale.edu/noonan_public/reilly2015/Track_Hub_me2/mm9/${timePoint}_k4me2_rep1.bw; mv ${timePoint}_k4me2_rep1.bw ${timePoint}_me2_rep1.bw
    wget http://noonan.ycga.yale.edu/noonan_public/reilly2015/Track_Hub_me2/mm9/${timePoint}_k4me2_rep2.bw; mv ${timePoint}_k4me2_rep2.bw ${timePoint}_me2_rep2.bw
done



############################## LIMB ######################################

# download bed files to /home/ak2267/project/EnhancerClasses/Limb
# H3K27ac
cd /home/ak2267/project/EnhancerClasses/Limb/hg19
for timePoint in E33 E41 E44 E47; do wget http://noonan.ycga.yale.edu/jcotney_www/Limb_hub/hg19/merge_${timePoint}_overlap.bed; done
cd /home/ak2267/project/EnhancerClasses/Limb/rheMac2
wget http://noonan.ycga.yale.edu/jcotney_www/Limb_hub/rheMac2/rm_e31_h3k27acEnrichedRegions0.00001.bed; mv rm_e31_h3k27acEnrichedRegions0.00001.bed preMerge_e31_h3k27acEnrichedRegions0.00001.bed
wget http://noonan.ycga.yale.edu/jcotney_www/Limb_hub/rheMac2/rm_e32_h3k27ac_1EnrichedRegions0.00001.bed; mv rm_e32_h3k27ac_1EnrichedRegions0.00001.bed preMerge_e32_h3k27acEnrichedRegions0.00001.bed
wget http://noonan.ycga.yale.edu/jcotney_www/Limb_hub/rheMac2/rm_e33_h3k27ac_1EnrichedRegions0.00001.bed; mv rm_e33_h3k27ac_1EnrichedRegions0.00001.bed preMerge_e33_h3k27acEnrichedRegions0.00001.bed
# get intersection of early time points in rhesus
bedtools intersect -a preMerge_e31_h3k27acEnrichedRegions0.00001.bed -b preMerge_e32_h3k27acEnrichedRegions0.00001.bed | bedtools intersect -a - -b preMerge_e33_h3k27acEnrichedRegions0.00001.bed > merge_e31_overlap.bed # only 20764 peaks - may be too stringent a strategy
wget http://noonan.ycga.yale.edu/jcotney_www/Limb_hub/rheMac2/rm_e36_h3k27ac_1EnrichedRegions0.00001.bed; mv rm_e36_h3k27ac_1EnrichedRegions0.00001.bed merge_e36_overlap.bed
cd /home/ak2267/project/EnhancerClasses/Limb/mm9
for timePoint in e10.5 e11.5 e12.5 e13.5; do wget http://noonan.ycga.yale.edu/jcotney_www/Limb_hub/mm9/merge_${timePoint}_overlap_0.00001.bed; mv merge_${timePoint}_overlap_0.00001.bed merge_${timePoint}_overlap.bed; done

# add peak number to 4th column of peak files to keep track of which ones lift
cd /home/ak2267/project/EnhancerClasses/Limb/
for species in hg19 rheMac2 mm9
do
    for file in ${species}/*.bed
	do name=`basename $file .bed`; awk '{ print $1"\t"$2"\t"$3"\tPeak_"NR}' ${file} > ${species}/${name}_named.bed
    done;done

# bed files are here:
/home/ak2267/project/EnhancerClasses/Limb/hg19/merge_${timePoint}_overlap_named.bed # timePoint = E33/E41/E44/E47
/home/ak2267/project/EnhancerClasses/Limb/rheMac2/merge_${timePoint}_overlap_named.bed # timePoint = e31/e36
/home/ak2267/project/EnhancerClasses/Limb/mm9/merge_${timePoint}_overlap_named.bed # timePoint = e10.5/e11.5/e12.5/e13.5

# get bigWig files
mkdir /gpfs/gibbs/pi/noonan/ak2267/NoonanLimb
cd /gpfs/gibbs/pi/noonan/ak2267/NoonanLimb
mkdir hg19; mkdir rheMac2; mkdir mm9

# hg19
cd /gpfs/gibbs/pi/noonan/ak2267/NoonanLimb/hg19
wget http://noonan.ycga.yale.edu/jcotney_www/Limb_hub/hg19/E33_N2102_h3k27ac_300_uniq.bw; mv E33_N2102_h3k27ac_300_uniq.bw E33_ac_rep1.bw
wget http://noonan.ycga.yale.edu/jcotney_www/Limb_hub/hg19/E33_N2104_h3k27ac_300_uniq.bw; mv E33_N2104_h3k27ac_300_uniq.bw E33_ac_rep2.bw
wget http://noonan.ycga.yale.edu/jcotney_www/Limb_hub/hg19/E41_N1699_h3k27ac_300_uniq.bw; mv E41_N1699_h3k27ac_300_uniq.bw E41_ac_rep1.bw
wget http://noonan.ycga.yale.edu/jcotney_www/Limb_hub/hg19/E41_11319_h3k27ac_300_uniq.bw; mv E41_11319_h3k27ac_300_uniq.bw E41_ac_rep2.bw
wget http://noonan.ycga.yale.edu/jcotney_www/Limb_hub/hg19/E44_11151_h3k27ac_300_uniq.bw; mv E44_11151_h3k27ac_300_uniq.bw E44_ac_rep1.bw
wget http://noonan.ycga.yale.edu/jcotney_www/Limb_hub/hg19/E44_11184_h3k27ac_1_300_uniq.bw; mv E44_11184_h3k27ac_1_300_uniq.bw E44_ac_rep2.bw
wget http://noonan.ycga.yale.edu/jcotney_www/Limb_hub/hg19/E47_11211_h3k27ac_300_uniq.bw; mv E47_11211_h3k27ac_300_uniq.bw E47_ac_rep1.bw
wget http://noonan.ycga.yale.edu/jcotney_www/Limb_hub/hg19/E47_11234_h3k27ac_1_300_uniq.bw; mv E47_11234_h3k27ac_1_300_uniq.bw E47_ac_rep2.bw

# rheMac2
cd /gpfs/gibbs/pi/noonan/ak2267/NoonanLimb/rheMac2
wget http://noonan.ycga.yale.edu/jcotney_www/Limb_hub/rheMac2/rm_e31_h3k27ac_300_uniq.bw; mv rm_e31_h3k27ac_300_uniq.bw e31_ac_rep1.bw
wget http://noonan.ycga.yale.edu/jcotney_www/Limb_hub/rheMac2/rm_e32_h3k27ac_1_300_uniq.bw; mv rm_e32_h3k27ac_1_300_uniq.bw e32_ac_rep1.bw
wget http://noonan.ycga.yale.edu/jcotney_www/Limb_hub/rheMac2/rm_e33_h3k27ac_1_300_uniq.bw; mv rm_e33_h3k27ac_1_300_uniq.bw e33_ac_rep1.bw
wget http://noonan.ycga.yale.edu/jcotney_www/Limb_hub/rheMac2/rm_e36_h3k27ac_1_300_uniq.bw; mv rm_e36_h3k27ac_1_300_uniq.bw e36_ac_rep1.bw

# mm9
cd /gpfs/gibbs/pi/noonan/ak2267/NoonanLimb/mm9
wget http://noonan.ycga.yale.edu/jcotney_www/Limb_hub/mm9/e10.5_1_h3k27ac_1_300_uniq.bw; mv e10.5_1_h3k27ac_1_300_uniq.bw e10.5_ac_rep1.bw
wget http://noonan.ycga.yale.edu/jcotney_www/Limb_hub/mm9/e10.5_2_h3k27ac_1_300_uniq.bw; mv e10.5_2_h3k27ac_1_300_uniq.bw e10.5_ac_rep2.bw
wget http://noonan.ycga.yale.edu/jcotney_www/Limb_hub/mm9/e10.5_3_h3k27ac_300_uniq.bw; mv e10.5_3_h3k27ac_300_uniq.bw e10.5_ac_rep3.bw
wget http://noonan.ycga.yale.edu/jcotney_www/Limb_hub/mm9/e11.5_1_h3k27ac_1_300_uniq.bw; mv e11.5_1_h3k27ac_1_300_uniq.bw e11.5_ac_rep1.bw
wget http://noonan.ycga.yale.edu/jcotney_www/Limb_hub/mm9/e11.5_2_h3k27ac_300_uniq.bw; mv e11.5_2_h3k27ac_300_uniq.bw e11.5_ac_rep2.bw
wget http://noonan.ycga.yale.edu/jcotney_www/Limb_hub/mm9/e12.5_1_h3k27ac_300_uniq.bw; mv e12.5_1_h3k27ac_300_uniq.bw e12.5_ac_rep1.bw
wget http://noonan.ycga.yale.edu/jcotney_www/Limb_hub/mm9/e12.5_2_h3k27ac_300_uniq.bw; mv e12.5_2_h3k27ac_300_uniq.bw e12.5_ac_rep2.bw
wget http://noonan.ycga.yale.edu/jcotney_www/Limb_hub/mm9/e13.5_1_h3k27ac_300_uniq.bw; mv e13.5_1_h3k27ac_300_uniq.bw e13.5_ac_rep1.bw
wget http://noonan.ycga.yale.edu/jcotney_www/Limb_hub/mm9/e13.5_2_h3k27ac_300_uniq.bw; mv e13.5_2_h3k27ac_300_uniq.bw e13.5_ac_rep2.bw


