# 2/10/23
# Choose example sites for rhesus vs mouse (rheMac10 vs mm39)

cd /gpfs/gibbs/pi/noonan/ak2267/speciesPairs/consensusCGIs/rheMac10_mm39

#########
# FIG 2 #
#########

# use this one for shared, 24 vs 25 CpGs
hg38 chr2	235630696	235631356	a_15628,b_1328 
rheMac10    chr12	123046553	123047211	a_15628,b_1328
mm39    chr1	89477819	89478434	a_15628,b_1328

# use this one for rhesus-only, 26 vs 8 CpGs
hg38 chr8	1453485	1454151	a_65786
rheMac10    chr8	1974029	1974695	a_65786
mm39 chr8	14704023	14704740	a_65786

## GET FASTA

# cat rheMac10_exampleSites.bed 
chr12	123046553	123047211	a_15628,b_1328 # 658 bp
chr8	1974029	1974695	a_65786 # 666 bp

# cat mm39_exampleSites.bed 
chr1	89477819	89478434	a_15628,b_1328 # 615 bp
chr8	14704023	14704740	a_65786 # 717 bp

bedtools getfasta -fi /gpfs/gibbs/pi/noonan/ak2267/genomes/rheMac10.fa -bed rheMac10_exampleSites.bed > rheMac10_exampleSites.fa
bedtools getfasta -fi /gpfs/gibbs/pi/noonan/ak2267/genomes/mm39.fa -bed mm39_exampleSites.bed > mm39_exampleSites.fa

## CALL CpG SITES
python findCpGs.py rheMac10_exampleSites.fa > rheMac10_exampleSites_CpGs.bed
python findCpGs.py mm39_exampleSites.fa > mm39_exampleSites_CpGs.bed

# view same region size so sites are comparable
# SHARED SITE
rheMac10  chr12 123046377   123047377
mm39    chr1    89477627    89478627

# RHESUS-ONLY SITE
rheMac10 chr8   1973864  1974864
mm39    chr8    14703882    14704882

## ADD AN ALIGNMENT POINT - ideally a centrally located CpG

# SHARED SITE alignment point
rheMac10    chr12   123046875   123046876
mm39    chr1    89478100    89478101

# RHESUS-ONLY SITE
rheMac10    chr8    1974418 1974419 
mm39    chr8    14704519    14704520

#########
# FIG 3 #
#########

hg38    chr13	36154825	36155556	a_31482
rheMac10  chr17	15063955	15064688	a_31482
mm39   chr3	55133889	55134555	a_31482

## GET FASTA

# cat rheMac10_Fig3_example.bed
chr17	15063955	15064688	a_31482 # 733 bp

# cat mm39_Fig3_example.bed 
chr3	55133889	55134555	a_31482 # 666 bp

bedtools getfasta -fi /gpfs/gibbs/pi/noonan/ak2267/genomes/rheMac10.fa -bed rheMac10_Fig3_example.bed > rheMac10_Fig3_example.fa
bedtools getfasta -fi /gpfs/gibbs/pi/noonan/ak2267/genomes/mm39.fa -bed mm39_Fig3_example.bed  > mm39_Fig3_example.fa

## CALL CpG SITES
python findCpGs.py rheMac10_Fig3_example.fa > rheMac10_Fig3_example_CpGs.bed # 40 CpGs
python findCpGs.py mm39_Fig3_example.fa > mm39_Fig3_example_CpGs.bed # 18 CpGs

# view same region size so sites are comparable
rheMac10    chr17   15062822    15065822
mm39    chr3    55132725    55135725

# temporarily change scale to 0 - 1 instead of 0 - 3 so peak fills the y axis

## ADD AN ALIGNMENT POINT
rheMac10    chr17   15064372    15064373
mm39 chr3   55134163    55134164
 
# HAVE TO FLIP MOUSE LOCUS!!!!




