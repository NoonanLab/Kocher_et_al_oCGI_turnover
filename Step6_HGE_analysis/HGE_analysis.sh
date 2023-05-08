# 1/28/23
# Purpose: analysis of oCGI species patterns in HGEs vs non-HGEs

# Improvements from previous pipeline:
# Add quantifications from bigWigs (using bigWigAverageOverBed) to enhancer table
# in order to more accurately control for enhancer strength in the background set of all enhancers
# Separate out time points
# Identify 3-way orthologous regions as another control - in this version just lifted one way to rhesus and mouse, no reciprocal lifts
# Add CpG numbers across entire peak rather than just CGI
# Add limb analysis, not just brain

# work here:
cd /home/ak2267/project/EnhancerClasses/HGE

# CGI files are here (in hg19, rheMac2, mm9):
/home/ak2267/genomes/CGI/UCSC_AL/${species}_CGIsAL.bed

# bed files are here (downloaded & modified in 220728_processNoonan.sh):
# brain
/home/ak2267/project/EnhancerClasses/${species}/${mark}/${file}_named.bed
# limb
/home/ak2267/project/EnhancerClasses/Limb/${species}/${file}_named.bed

# make directories for intermediate files
for tissue in brain limb
do
    mkdir ${tissue}
    for species in hg19 rheMac2 mm9
    do
        mkdir ${tissue}/${species}
        mkdir ${tissue}/${species}/ac
        mkdir ${tissue}/${species}/me2 # remove from limb directory
    done
done

### intersect peaks with CGIs in each species
cd /home/ak2267/project/EnhancerClasses/HGE

for species in hg19 rheMac2 mm9
do
    # brain (ac and me2)
    for mark in ac me2
    do
	    for file in /home/ak2267/project/EnhancerClasses/${species}/${mark}/*_named.bed
	    do name=`basename $file _named.bed`; bedtools intersect -wao -a ${file} -b /home/ak2267/genomes/CGI/UCSC_AL/${species}_CGIsAL.bed > brain/${species}/${mark}/${name}_intersectCGIs.txt
	    done;done
	# limb (ac only)
	mark=ac
	for file in /home/ak2267/project/EnhancerClasses/Limb/${species}/*_named.bed
	do name=`basename $file _named.bed`; bedtools intersect -wao -a ${file} -b /home/ak2267/genomes/CGI/UCSC_AL/${species}_CGIsAL.bed > limb/${species}/${mark}/${name}_intersectCGIs.txt
    done
done

### intersect human merged peaks with HGEs (enhancers) and HGPs (promoters)
# download gain files for brain ac
cd /home/ak2267/project/EnhancerClasses/hg19/Gains/ac
for timePoint in CS16 CS23 F2F F2O; do wget http://noonan.ycga.yale.edu/noonan_public/reilly2015/Track_Hub_ac/hg19/${timePoint}_overlap_hs_gain_enhancer.bb; done
for timePoint in CS16 CS23 F2F F2O; do wget http://noonan.ycga.yale.edu/noonan_public/reilly2015/Track_Hub_ac/hg19/${timePoint}_overlap_hs_gain_promoter.bb; done
# convert bigBed to bed and remove bigBeds
for file in $(ls *.bb); do name=`basename $file .bb`; bigBedToBed ${file} ${name}.bed ; done
rm *.bb
# download gain files for brain me2
cd /home/ak2267/project/EnhancerClasses/hg19/Gains/me2
for timePoint in CS16 CS23 F2F F2O; do wget http://noonan.ycga.yale.edu/noonan_public/reilly2015/Track_Hub_me2/hg19/${timePoint}_me2_overlap_hs_gain_enhancer.bb; done
for timePoint in CS16 CS23 F2F F2O; do wget http://noonan.ycga.yale.edu/noonan_public/reilly2015/Track_Hub_me2/hg19/${timePoint}_me2_overlap_hs_gain_promoter.bb; done
# convert bigBed to bed and remove bigBeds
for file in $(ls *.bb); do name=`basename $file .bb`; bigBedToBed ${file} ${name}.bed ; done
rm *.bb
# download gain files for limb ac
cd /home/ak2267/project/EnhancerClasses/hg19/Gains/limbAc
for timePoint in E33 E41 E44 E47; do wget http://noonan.ycga.yale.edu/jcotney_www/Limb_hub/hg19/bhP_0.001_log2FC_0.5849625_merge_all_hs_gain_${timePoint}.bed; done
for timePoint in E33 E41 E44 E47; do mv bhP_0.001_log2FC_0.5849625_merge_all_hs_gain_${timePoint}.bed merge_all_hs_gain_${timePoint}.bed; done

# intersect all peaks with Gains
cd /home/ak2267/project/EnhancerClasses/HGE
# brain
for timePoint in CS16 CS23 F2F F2O
do
    bedtools intersect -wao -a /home/ak2267/project/EnhancerClasses/hg19/ac/merge_${timePoint}_overlap_named.bed -b /home/ak2267/project/EnhancerClasses/hg19/Gains/ac/${timePoint}_overlap_hs_gain_enhancer.bed > brain/hg19/ac/merge_${timePoint}_overlap_intersectHGEs.txt
    bedtools intersect -wao -a /home/ak2267/project/EnhancerClasses/hg19/ac/merge_${timePoint}_overlap_named.bed -b /home/ak2267/project/EnhancerClasses/hg19/Gains/ac/${timePoint}_overlap_hs_gain_promoter.bed > brain/hg19/ac/merge_${timePoint}_overlap_intersectHGPs.txt
    bedtools intersect -wao -a /home/ak2267/project/EnhancerClasses/hg19/me2/merge_${timePoint}_me2_overlap_named.bed -b /home/ak2267/project/EnhancerClasses/hg19/Gains/me2/${timePoint}_me2_overlap_hs_gain_enhancer.bed > brain/hg19/me2/merge_${timePoint}_overlap_intersectHGEs.txt
    bedtools intersect -wao -a /home/ak2267/project/EnhancerClasses/hg19/me2/merge_${timePoint}_me2_overlap_named.bed -b /home/ak2267/project/EnhancerClasses/hg19/Gains/me2/${timePoint}_me2_overlap_hs_gain_promoter.bed > brain/hg19/me2/merge_${timePoint}_overlap_intersectHGPs.txt
done
# limb
for timePoint in E33 E41 E44 E47
do
    bedtools intersect -wao -a /home/ak2267/project/EnhancerClasses/Limb/hg19/merge_${timePoint}_overlap_named.bed -b /home/ak2267/project/EnhancerClasses/hg19/Gains/limbAc/merge_all_hs_gain_${timePoint}.bed > limb/hg19/ac/merge_${timePoint}_overlap_intersectGains.txt
done


#### remove those that I annotate as a promoter or other feature
# lift my gene annotations from hg38 to hg19 - done in 220511_HGEs.sh
cd /home/ak2267/project/EnhancerClasses/featureAnnotations
for file in $(ls /home/ak2267/genomes/RefSeq/featureAnnotations/hg38*); do name=`basename $file .bed`; liftOver ${file} /home/ak2267/genomes/chain/hg38ToHg19.over.chain.gz ${name}_LOhg19.bed unMapped; done

# intersect human merged peaks with protein-coding promoters
cd /home/ak2267/project/EnhancerClasses/HGE
# brain
for timePoint in CS16 CS23 F2F F2O
do
    bedtools intersect -wa -u -a /home/ak2267/project/EnhancerClasses/hg19/ac/merge_${timePoint}_overlap_named.bed \
	     -b ../featureAnnotations/hg38_1_proteinCoding_promoters_merged_LOhg19.bed > brain/hg19/ac/merge_${timePoint}_overlap_intersectProteinCodingPromoters.bed
    bedtools intersect -wa -u -a /home/ak2267/project/EnhancerClasses/hg19/me2/merge_${timePoint}_me2_overlap_named.bed \
	     -b ../featureAnnotations/hg38_1_proteinCoding_promoters_merged_LOhg19.bed > brain/hg19/me2/merge_${timePoint}_overlap_intersectProteinCodingPromoters.bed
done
# limb
for timePoint in E33 E41 E44 E47
do
    bedtools intersect -wa -u -a /home/ak2267/project/EnhancerClasses/Limb/hg19/merge_${timePoint}_overlap_named.bed \
        -b ../featureAnnotations/hg38_1_proteinCoding_promoters_merged_LOhg19.bed > limb/hg19/ac/merge_${timePoint}_overlap_intersectProteinCodingPromoters.bed
done


## intersect human merged peaks with enhancer space (intronic & intergenic)
# brain
for timePoint in CS16 CS23 F2F F2O
do
    bedtools intersect -v -a /home/ak2267/project/EnhancerClasses/hg19/ac/merge_${timePoint}_overlap_named.bed \
	 -b ../featureAnnotations/hg38_1_proteinCoding_promoters_merged_LOhg19.bed \
	 -b ../featureAnnotations/hg38_2_proteinCoding_exons_merged_LOhg19.bed \
	 -b ../featureAnnotations/hg38_3_lncRNA_merged_LOhg19.bed \
	 -b ../featureAnnotations/hg38_4_ncRNA_merged_LOhg19.bed \
	 -b ../featureAnnotations/hg38_5_pseudogene_merged_LOhg19.bed \
	 -b ../featureAnnotations/hg38_6_unknown_merged_LOhg19.bed > brain/hg19/ac/merge_${timePoint}_overlap_intersectIntronicIntergenic.bed
    bedtools intersect -v -a /home/ak2267/project/EnhancerClasses/hg19/me2/merge_${timePoint}_me2_overlap_named.bed \
	 -b ../featureAnnotations/hg38_1_proteinCoding_promoters_merged_LOhg19.bed \
	 -b ../featureAnnotations/hg38_2_proteinCoding_exons_merged_LOhg19.bed \
	 -b ../featureAnnotations/hg38_3_lncRNA_merged_LOhg19.bed \
	 -b ../featureAnnotations/hg38_4_ncRNA_merged_LOhg19.bed \
	 -b ../featureAnnotations/hg38_5_pseudogene_merged_LOhg19.bed \
	 -b ../featureAnnotations/hg38_6_unknown_merged_LOhg19.bed > brain/hg19/me2/merge_${timePoint}_overlap_intersectIntronicIntergenic.bed
done
# limb
for timePoint in E33 E41 E44 E47
do
    bedtools intersect -v -a /home/ak2267/project/EnhancerClasses/Limb/hg19/merge_${timePoint}_overlap_named.bed \
	 -b ../featureAnnotations/hg38_1_proteinCoding_promoters_merged_LOhg19.bed \
	 -b ../featureAnnotations/hg38_2_proteinCoding_exons_merged_LOhg19.bed \
	 -b ../featureAnnotations/hg38_3_lncRNA_merged_LOhg19.bed \
	 -b ../featureAnnotations/hg38_4_ncRNA_merged_LOhg19.bed \
	 -b ../featureAnnotations/hg38_5_pseudogene_merged_LOhg19.bed \
	 -b ../featureAnnotations/hg38_6_unknown_merged_LOhg19.bed > limb/hg19/ac/merge_${timePoint}_overlap_intersectIntronicIntergenic.bed
done

# lift over human peaks to rhesus and to mouse
# use -minMatch=0.1 because that captures all human gains (can't find or reverse-engineer the setting used in Reilly, Yin 2015)
cd /home/ak2267/project/EnhancerClasses/HGE
# brain
for timePoint in CS16 CS23 F2F F2O
do
    liftOver -minMatch=0.1 /home/ak2267/project/EnhancerClasses/hg19/ac/merge_${timePoint}_overlap_named.bed /home/ak2267/genomes/chain/hg19ToRheMac2.over.chain.gz brain/hg19/ac/merge_${timePoint}_overlap_LOrheMac2.bed unMapped
    liftOver -minMatch=0.1 /home/ak2267/project/EnhancerClasses/hg19/me2/merge_${timePoint}_me2_overlap_named.bed /home/ak2267/genomes/chain/hg19ToRheMac2.over.chain.gz brain/hg19/me2/merge_${timePoint}_overlap_LOrheMac2.bed unMapped
    liftOver -minMatch=0.1 /home/ak2267/project/EnhancerClasses/hg19/ac/merge_${timePoint}_overlap_named.bed /home/ak2267/genomes/chain/hg19ToMm9.over.chain.gz brain/hg19/ac/merge_${timePoint}_overlap_LOmm9.bed unMapped
    liftOver -minMatch=0.1 /home/ak2267/project/EnhancerClasses/hg19/me2/merge_${timePoint}_me2_overlap_named.bed /home/ak2267/genomes/chain/hg19ToMm9.over.chain.gz brain/hg19/me2/merge_${timePoint}_overlap_LOmm9.bed unMapped
done
# limb
for timePoint in E33 E41 E44 E47
do
    liftOver -minMatch=0.1 /home/ak2267/project/EnhancerClasses/Limb/hg19/merge_${timePoint}_overlap_named.bed /home/ak2267/genomes/chain/hg19ToRheMac2.over.chain.gz limb/hg19/ac/merge_${timePoint}_overlap_LOrheMac2.bed unMapped
    liftOver -minMatch=0.1 /home/ak2267/project/EnhancerClasses/Limb/hg19/merge_${timePoint}_overlap_named.bed /home/ak2267/genomes/chain/hg19ToMm9.over.chain.gz limb/hg19/ac/merge_${timePoint}_overlap_LOmm9.bed unMapped
done

# lift back to human
# brain
for timePoint in CS16 CS23 F2F F2O
do
    liftOver -minMatch=0.1 brain/hg19/ac/merge_${timePoint}_overlap_LOrheMac2.bed /home/ak2267/genomes/chain/rheMac2ToHg19.over.chain.gz brain/hg19/ac/merge_${timePoint}_overlap_LOrheMac2_backToHg19.bed unMapped
    liftOver -minMatch=0.1 brain/hg19/me2/merge_${timePoint}_overlap_LOrheMac2.bed /home/ak2267/genomes/chain/rheMac2ToHg19.over.chain.gz brain/hg19/me2/merge_${timePoint}_overlap_LOrheMac2_backToHg19.bed unMapped
    liftOver -minMatch=0.1 brain/hg19/ac/merge_${timePoint}_overlap_LOmm9.bed /home/ak2267/genomes/chain/mm9ToHg19.over.chain.gz brain/hg19/ac/merge_${timePoint}_overlap_LOmm9_backToHg19.bed unMapped
    liftOver -minMatch=0.1 brain/hg19/me2/merge_${timePoint}_overlap_LOmm9.bed /home/ak2267/genomes/chain/mm9ToHg19.over.chain.gz brain/hg19/me2/merge_${timePoint}_overlap_LOmm9_backToHg19.bed unMapped
done
# limb
for timePoint in E33 E41 E44 E47
do
    liftOver -minMatch=0.1 limb/hg19/ac/merge_${timePoint}_overlap_LOrheMac2.bed /home/ak2267/genomes/chain/rheMac2ToHg19.over.chain.gz limb/hg19/ac/merge_${timePoint}_overlap_LOrheMac2_backToHg19.bed unMapped
    liftOver -minMatch=0.1 limb/hg19/ac/merge_${timePoint}_overlap_LOmm9.bed /home/ak2267/genomes/chain/mm9ToHg19.over.chain.gz limb/hg19/ac/merge_${timePoint}_overlap_LOmm9_backToHg19.bed unMapped
done

# restrict to reciprocal lifts
# brain
for timePoint in CS16 CS23 F2F F2O
do
    echo 'starting'
    for mark in ac me2
    do
        for species in rheMac2 mm9
        do
            bedtools intersect -wao -a brain/hg19/${mark}/merge_${timePoint}_overlap_LO${species}_backToHg19.bed -b /home/ak2267/project/EnhancerClasses/hg19/${mark}/merge_${timePoint}_*overlap_named.bed > brain/hg19/${mark}/${species}_${timePoint}_checkMapping.txt
            python /gpfs/gibbs/pi/noonan/ak2267/speciesPairs/consensusCGIs/Scripts/restrictToSitesThatMapBack.py brain/hg19/${mark}/${species}_${timePoint}_checkMapping.txt > brain/hg19/${mark}/hg19_${mark}_${timePoint}_sitesThatMapTo${species^}_inHg19.bed
            python /gpfs/gibbs/pi/noonan/ak2267/speciesPairs/consensusCGIs/Scripts/restrictToLO.py brain/hg19/${mark}/merge_${timePoint}_overlap_LO${species}.bed brain/hg19/${mark}/hg19_${mark}_${timePoint}_sitesThatMapTo${species^}_inHg19.bed > brain/hg19/${mark}/hg19_${mark}_${timePoint}_sitesThatMapTo${species^}_in${species^}.bed
        done
    done
done
# limb
for timePoint in E33 E41 E44 E47
do
    mark=ac
    for species in rheMac2 mm9
    do
        bedtools intersect -wao -a limb/hg19/${mark}/merge_${timePoint}_overlap_LO${species}_backToHg19.bed -b /home/ak2267/project/EnhancerClasses/Limb/hg19/merge_${timePoint}_*overlap_named.bed > limb/hg19/${mark}/${species}_${timePoint}_checkMapping.txt
        python /gpfs/gibbs/pi/noonan/ak2267/speciesPairs/consensusCGIs/Scripts/restrictToSitesThatMapBack.py limb/hg19/${mark}/${species}_${timePoint}_checkMapping.txt > limb/hg19/${mark}/hg19_${mark}_${timePoint}_sitesThatMapTo${species^}_inHg19.bed
        python /gpfs/gibbs/pi/noonan/ak2267/speciesPairs/consensusCGIs/Scripts/restrictToLO.py limb/hg19/${mark}/merge_${timePoint}_overlap_LO${species}.bed limb/hg19/${mark}/hg19_${mark}_${timePoint}_sitesThatMapTo${species^}_inHg19.bed > limb/hg19/${mark}/hg19_${mark}_${timePoint}_sitesThatMapTo${species^}_in${species^}.bed
    done    
done


# intersect with CGIs in rhesus and mouse
cd /home/ak2267/project/EnhancerClasses/HGE
for species in rheMac2 mm9
do
    # brain
    for timePoint in CS16 CS23 F2F F2O
    do 
	for mark in ac me2
	do
	    bedtools intersect -wao -a brain/hg19/${mark}/hg19_${mark}_${timePoint}_sitesThatMapTo${species^}_in${species^}.bed \
	    -b /home/ak2267/genomes/CGI/UCSC_AL/${species}_CGIsAL.bed > brain/hg19/${mark}/hg19_${mark}_${timePoint}_sitesThatMapTo${species^}_intersect${species^}CGIs.bed
	done;done
	# limb
	for timePoint in E33 E41 E44 E47
	do
	mark=ac
	bedtools intersect -wao -a limb/hg19/${mark}/hg19_${mark}_${timePoint}_sitesThatMapTo${species^}_in${species^}.bed \
	    -b /home/ak2267/genomes/CGI/UCSC_AL/${species}_CGIsAL.bed > limb/hg19/${mark}/hg19_${mark}_${timePoint}_sitesThatMapTo${species^}_intersect${species^}CGIs.bed
	done;done
	

# bigWigs are here:
# /gpfs/gibbs/pi/noonan/ak2267/NoonanBrain/hg19/ac
# /gpfs/gibbs/pi/noonan/ak2267/NoonanBrain/hg19/me2
# /gpfs/gibbs/pi/noonan/ak2267/NoonanLimb/hg19/

# bed files for peaks (merged across replicates) were downloaded to /home/ak2267/project/EnhancerClasses/hg19 above
# /home/ak2267/project/EnhancerClasses/hg19/ac/merge_${timePoint}_overlap_named.bed
# /home/ak2267/project/EnhancerClasses/hg19/me2/merge_${timePoint}_me2_overlap_named.bed
# /home/ak2267/project/EnhancerClasses/Limb/hg19/merge_${timePoint}_overlap_named.bed

## QUANTIFY signal in peaks in each replicate using bigWigAverageOverBed
# run bigWigAverageOverBed
# The output columns are:
   #name - name field from bed, which should be unique
   #size - size of bed (sum of exon sizes)
   #covered - # bases within exons covered by bigWig
   #sum - sum of values over all bases covered
   #mean0 - average over bases with non-covered bases counting as zeroes
   #mean - average over just covered bases
cd /home/ak2267/project/EnhancerClasses/quant/brain/hg19/
# CS16
bigWigAverageOverBed /home/ak2267/scratch60/HGE/CS16/GSM1554660_Hu_7pcw_H3K27ac_rep1.bw \
		     /home/ak2267/project/EnhancerClasses/hg19/ac/merge_CS16_overlap_named.bed ac/CS16_ac_quant_rep1.tab
bigWigAverageOverBed /home/ak2267/scratch60/HGE/CS16/GSM1554661_Hu_7pcw_H3K27ac_rep2.bw \
		     /home/ak2267/project/EnhancerClasses/hg19/ac/merge_CS16_overlap_named.bed ac/CS16_ac_quant_rep2.tab
bigWigAverageOverBed /home/ak2267/scratch60/HGE/CS16/GSM1554662_Hu_7pcw_H3K4me2_rep1.bw \
		     /home/ak2267/project/EnhancerClasses/hg19/me2/merge_CS16_me2_overlap_named.bed me2/CS16_me2_quant_rep1.tab
bigWigAverageOverBed /home/ak2267/scratch60/HGE/CS16/GSM1554663_Hu_7pcw_H3K4me2_rep2.bw \
		     /home/ak2267/project/EnhancerClasses/hg19/me2/merge_CS16_me2_overlap_named.bed me2/CS16_me2_quant_rep2.tab
# CS23
bigWigAverageOverBed /home/ak2267/scratch60/HGE/CS23/GSM1554666_Hu_8_5pcw_H3K27ac_rep1.bw \
		     /home/ak2267/project/EnhancerClasses/hg19/ac/merge_CS23_overlap_named.bed ac/CS23_ac_quant_rep1.tab
bigWigAverageOverBed /home/ak2267/scratch60/HGE/CS23/GSM1554667_Hu_8_5pcw_H3K27ac_rep2.bw \
		     /home/ak2267/project/EnhancerClasses/hg19/ac/merge_CS23_overlap_named.bed ac/CS23_ac_quant_rep2.tab
bigWigAverageOverBed /home/ak2267/scratch60/HGE/CS23/GSM1554668_Hu_8_5pcw_H3K4me2_rep1.bw \
		     /home/ak2267/project/EnhancerClasses/hg19/me2/merge_CS23_me2_overlap_named.bed me2/CS23_me2_quant_rep1.tab
bigWigAverageOverBed /home/ak2267/scratch60/HGE/CS23/GSM1554669_Hu_8_5pcw_H3K4me2_rep2.bw \
		     /home/ak2267/project/EnhancerClasses/hg19/me2/merge_CS23_me2_overlap_named.bed me2/CS23_me2_quant_rep2.tab

# F2F
bigWigAverageOverBed /home/ak2267/scratch60/HGE/F2F/GSM1554672_Hu_12Fpcw_H3K27ac_rep1.bw \
		     /home/ak2267/project/EnhancerClasses/hg19/ac/merge_F2F_overlap_named.bed ac/F2F_ac_quant_rep1.tab
bigWigAverageOverBed /home/ak2267/scratch60/HGE/F2F/GSM1554673_Hu_12Fpcw_H3K27ac_rep2.bw \
		     /home/ak2267/project/EnhancerClasses/hg19/ac/merge_F2F_overlap_named.bed ac/F2F_ac_quant_rep2.tab
bigWigAverageOverBed /home/ak2267/scratch60/HGE/F2F/GSM1554674_Hu_12Fpcw_H3K4me2_rep1.bw \
		     /home/ak2267/project/EnhancerClasses/hg19/me2/merge_F2F_me2_overlap_named.bed me2/F2F_me2_quant_rep1.tab
bigWigAverageOverBed /home/ak2267/scratch60/HGE/F2F/GSM1554675_Hu_12Fpcw_H3K4me2_rep2.bw \
		     /home/ak2267/project/EnhancerClasses/hg19/me2/merge_F2F_me2_overlap_named.bed me2/F2F_me2_quant_rep2.tab

# F2O
bigWigAverageOverBed /home/ak2267/scratch60/HGE/F2O/GSM1554678_Hu_12Opcw_H3K27ac_rep1.bw \
    /home/ak2267/project/EnhancerClasses/hg19/ac/merge_F2O_overlap_named.bed ac/F2O_ac_quant_rep1.tab
bigWigAverageOverBed /home/ak2267/scratch60/HGE/F2O/GSM1554679_Hu_12Opcw_H3K27ac_rep2.bw \
	/home/ak2267/project/EnhancerClasses/hg19/ac/merge_F2O_overlap_named.bed ac/F2O_ac_quant_rep2.tab
bigWigAverageOverBed /home/ak2267/scratch60/HGE/F2O/GSM1554680_Hu_12Opcw_H3K4me2_rep1.bw \
	/home/ak2267/project/EnhancerClasses/hg19/me2/merge_F2O_me2_overlap_named.bed me2/F2O_me2_quant_rep1.tab
bigWigAverageOverBed /home/ak2267/scratch60/HGE/F2O/GSM1554681_Hu_12Opcw_H3K4me2_rep2.bw \
	/home/ak2267/project/EnhancerClasses/hg19/me2/merge_F2O_me2_overlap_named.bed me2/F2O_me2_quant_rep2.tab

# LIMB
cd /home/ak2267/project/EnhancerClasses/quant/limb/hg19/
# E33
for timePoint in E33 E41 E44 E47
do
    for rep in rep1 rep2
    do
        echo $timePoint $rep
        bigWigAverageOverBed /gpfs/gibbs/pi/noonan/ak2267/NoonanLimb/hg19/${timePoint}_ac_${rep}.bw \
        /home/ak2267/project/EnhancerClasses/Limb/hg19/merge_${timePoint}_overlap_named.bed ac/${timePoint}_ac_quant_${rep}.tab
    done
done

## Count CpGs in peaks, not just CGIs, using faCount
# bedtools getfasta [OPTIONS] -fi <fasta> -bed <bed/gff/vcf>
# genomes (unmasked) are already downloaded here:
# /gpfs/gibbs/pi/noonan/ak2267/genomes/${species}.fa

cd /home/ak2267/project/EnhancerClasses/HGE
mkdir CpG; cd /home/ak2267/project/EnhancerClasses/HGE/CpG
mkdir brain; mkdir limb
for species in hg19 rheMac2 mm9
do
    mkdir brain/${species}; mkdir brain/${species}/ac; mkdir brain/${species}/me2
    mkdir limb/${species}; mkdir limb/${species}/ac
done

# brain
for timePoint in CS16 CS23 F2F F2O
do
    echo $timePoint
    bedtools getfasta -nameOnly -fi /gpfs/gibbs/pi/noonan/ak2267/genomes/hg19.fa \
    -bed /home/ak2267/project/EnhancerClasses/hg19/ac/merge_${timePoint}_overlap_named.bed > brain/hg19/ac/${timePoint}.fa
    bedtools getfasta -nameOnly -fi /gpfs/gibbs/pi/noonan/ak2267/genomes/hg19.fa \
    -bed /home/ak2267/project/EnhancerClasses/hg19/me2/merge_${timePoint}_me2_overlap_named.bed > brain/hg19/me2/${timePoint}.fa
done
# limb
for timePoint in E33 E41 E44 E47
do
    echo $timePoint
    bedtools getfasta -nameOnly -fi /gpfs/gibbs/pi/noonan/ak2267/genomes/hg19.fa \
	    -bed /home/ak2267/project/EnhancerClasses/Limb/hg19/merge_${timePoint}_overlap_named.bed > limb/hg19/ac/${timePoint}.fa
done
														     
# faCount file(s).fa
# brain
for mark in ac me2; do for timePoint in CS16 CS23 F2F F2O; do echo $timePoint; faCount brain/hg19/${mark}/${timePoint}.fa > brain/hg19/${mark}/${timePoint}.faCount; done; done
# limb
for mark in ac; do for timePoint in E33 E41 E44 E47; do echo $timePoint; faCount limb/hg19/${mark}/${timePoint}.fa > limb/hg19/${mark}/${timePoint}.faCount; done; done

# do the same for RHESUS and MOUSE
cd /home/ak2267/project/EnhancerClasses/HGE/CpG
# brain
for mark in ac me2
do
    for timePoint in CS16 CS23 F2F F2O
    do
	for species in rheMac2 mm9
	do
	    echo $mark $timePoint $species
	    bedtools getfasta -nameOnly -fi /gpfs/gibbs/pi/noonan/ak2267/genomes/${species}.fa \
		     -bed /home/ak2267/project/EnhancerClasses/HGE/brain/hg19/${mark}/hg19_${mark}_${timePoint}_sitesThatMapTo${species^}_in${species^}.bed > brain/${species}/${mark}/${timePoint}.fa
	    faCount brain/${species}/${mark}/${timePoint}.fa > brain/${species}/${mark}/${timePoint}.faCount
	done;done;done
# limb
mark=ac
for timePoint in E33 E41 E44 E47
do
	for species in rheMac2 mm9
	do 
        echo $mark $timePoint $species
	    bedtools getfasta -nameOnly -fi /gpfs/gibbs/pi/noonan/ak2267/genomes/${species}.fa \
		     -bed /home/ak2267/project/EnhancerClasses/HGE/limb/hg19/${mark}/hg19_${mark}_${timePoint}_sitesThatMapTo${species^}_in${species^}.bed > limb/${species}/${mark}/${timePoint}.fa
	    faCount limb/${species}/${mark}/${timePoint}.fa > limb/${species}/${mark}/${timePoint}.faCount
	done;done

# run modified python script to incorporate bigWig quantification and CpG counts
cd /home/ak2267/project/EnhancerClasses/HGE
# brain
for mark in ac me2; do for timePoint in CS16 CS23 F2F F2O; do echo $mark $timePoint; python makeHGEtable_quant.py brain ${mark} ${timePoint} > HGEtable_brain_${mark}_${timePoint}.txt; done; done
# limb
for mark in ac; do for timePoint in E33 E41 E44 E47; do echo $mark $timePoint; python makeHGEtable_quant.py limb ${mark} ${timePoint} > HGEtable_limb_${mark}_${timePoint}.txt; done; done
							      
# download
mkdir HGEtable
mv HGEtable*.txt HGEtable
tar -zcvf HGEtable.gz HGEtable

cd /Users/acadiak/Desktop/CGI/EnhancerClasses/R_HGEs
scp ak2267@ruddle.hpc.yale.edu:/home/ak2267/project/EnhancerClasses/HGE/HGEtable.gz .


############ RUN RESAMPLING TESTS ON THE CLUSTER ############

cd /home/ak2267/project/EnhancerClasses/HGE/permutation

tissue=brain
for mark in ac me2
do
    for timePoint in CS16 CS23 F2F F2O
    do
        echo 'cd /home/ak2267/project/EnhancerClasses/HGE/permutation ; source ~/.bashrc ; source ~/.bash_profile ; module load R ; Rscript HGEs_resampling.R '${tissue}'_'${mark}'_'${timePoint}
    done; done >> 230215_resampling_HGE_jobFile.txt
tissue=limb
mark=ac
for timePoint in E33 E41 E44 E47
do
    echo 'cd /home/ak2267/project/EnhancerClasses/HGE/permutation ; source ~/.bashrc ; source ~/.bash_profile ; module load R ; Rscript HGEs_resampling.R '${tissue}'_'${mark}'_'${timePoint}
    done >> 230215_resampling_HGE_jobFile.txt
    
dsq --job-file 230215_resampling_HGE_jobFile.txt --mem-per-cpu 5G -c 1 --mail-type FAIL,END
sbatch dsq-230215_resampling_HGE_jobFile-2023-02-15.sh # 21178156
sbatch dsq-230215_resampling_HGE_jobFile-2023-03-09.sh # 21677366

# a few cases have slightly fewer than 20,000 resamples

# Summarize results of resampling in tables for use in R

for fileInfo in brain_ac_CS16 brain_ac_CS23 brain_ac_F2F brain_ac_F2O brain_me2_CS16 brain_me2_CS23 brain_me2_F2F brain_me2_F2O limb_ac_E33 limb_ac_E41 limb_ac_E44 limb_ac_E47
do
    python summarizeResamplingHGEs.py ${fileInfo}
done >> HGEs_resamplingSummary.txt

# download
cd /Users/acadiak/Desktop/CGI/EnhancerClasses/R_HGEs
scp ak2267@ruddle.hpc.yale.edu:/home/ak2267/project/EnhancerClasses/HGE/permutation/HGEs_resamplingSummary.txt .

# download intermediate files for CS23 brain for use in making supplemental figure
cd /Users/acadiak/Desktop/CGI/EnhancerClasses/R_HGEs
scp ak2267@ruddle.hpc.yale.edu:/home/ak2267/project/EnhancerClasses/HGE/permutation/brain_ac_CS23_observed_HGEs_CGI_patterns.txt .
scp ak2267@ruddle.hpc.yale.edu:/home/ak2267/project/EnhancerClasses/HGE/permutation/brain_ac_CS23_resampled_HGEs_CGI_patterns.txt .



