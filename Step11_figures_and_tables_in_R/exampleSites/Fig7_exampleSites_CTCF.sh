# 2/24/23
# Purpose: choose and visualize example sites for rhesus vs mouse CTCF binding

cd /gpfs/gibbs/pi/noonan/ak2267/speciesPairs/consensusCGIs/rheMac10_mm39

#########
# FIG 7 #
#########

# peaks are here
/gpfs/gibbs/pi/noonan/ak2267/LiverTF/peaks/CTCF/intersection/${species}_liver_CTCF_intersection.bed

# choose sites using code at bottom of Fig3-4_grids.R

# use this one for shared, 60 vs 52 CpGs
hg38 chr10	100670822	100671675	a_70865,b_28241
rheMac10    chr9	101143356	101144206	a_70865,b_28241
mm39 chr19	44672403	44673158	a_70865,b_28241

# use this one for rhesus-only, 25 vs 6 CpGs
hg38 chr11	116357956	116358338	a_21974
rheMac10    chr14	109432692	109433073	a_21974
mm39    chr9	46574650	46575026	a_21974


## GET FASTA

# cat rheMac10_exampleSites_CTCF.bed 
chr9	101143356	101144206	a_70865,b_28241 # 850 bp
chr14	109432692	109433073	a_21974 # 381 bp

# cat mm39_exampleSites_CTCF.bed 
chr19	44672403	44673158	a_70865,b_28241 # 755 bp
chr9	46574650	46575026	a_21974 # 376 bp

bedtools getfasta -fi /gpfs/gibbs/pi/noonan/ak2267/genomes/rheMac10.fa -bed rheMac10_exampleSites_CTCF.bed > rheMac10_exampleSites_CTCF.fa
bedtools getfasta -fi /gpfs/gibbs/pi/noonan/ak2267/genomes/mm39.fa -bed mm39_exampleSites_CTCF.bed > mm39_exampleSites_CTCF.fa

## CALL CpG SITES
python findCpGs.py rheMac10_exampleSites_CTCF.fa > rheMac10_exampleSites_CTCF_CpGs.bed
python findCpGs.py mm39_exampleSites_CTCF.fa > mm39_exampleSites_CTCF_CpGs.bed

# view same region size so sites are comparable (2kb)
# SHARED SITE
rheMac10  chr9  101142766   101144766
mm39    chr19   44671773    44673773

# RHESUS-ONLY SITE
rheMac10 chr14  109431878   109433878
mm39    chr9    46573888    46575887

## ADD AN ALIGNMENT POINT - ideally a centrally located CpG
# this is a 1-bp bed interval that I will align visually in Illustrator

# SHARED SITE alignment point - loci go in same direction
rheMac10    chr9    101143721   101143722
mm39    chr19   44672789    44672790

# RHESUS-ONLY SITE - loci go in OPPOSITE DIRECTIONS SO FLIP MOUSE LOCUS
rheMac10    chr14   109432867   109432868 
mm39    chr9    46574853    46574854
