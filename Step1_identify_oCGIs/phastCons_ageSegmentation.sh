#############
# phastCons #
#############

# 7/29/22
# Purpose: download 100-way phastCons elements in human genome for use in downstream pipeline with oCGI analysis

# download human 100-way phastCons elements from the Table Browser
# hg38 --> Comparative Genomics --> Conservation --> 100 Vert. El (phastConsElements100way)
# download as phastConsElements100way_hg38.txt
# put into ~/genomes/phastCons on Ruddle & cut into bed files
# updated on 8/9/22 to use LOD score in column 4 instead of other score
# https://genome.ucsc.edu/cgi-bin/hgTrackUi?db=hg19&c=chrX&g=cons46way

cd /home/ak2267/genomes/phastCons
cut -f 2,3,4,5 ~/genomes/phastCons/phastConsElements100way_hg38.txt | tail -10350729 > ~/genomes/phastCons/phastConsElements100way_hg38.bed

liftOver phastConsElements100way_hg38.bed /home/ak2267/genomes/chain/hg38ToHg19.over.chain.gz phastConsElements100way_hg19.bed unMapped


####################
# age segmentation #
####################

# 8/9/22
# Download age segmentation map from Emera, Yin et al 2016 for use in pipeline

mkdir /home/ak2267/project/EnhancerAge

# Download segmentation map from Emera & Yin 2016 and convert to bed file
cd /home/ak2267/project/EnhancerAge
wget http://noonan.ycga.yale.edu/noonan_public/emera2016/Age/hg19/hg19_nh.chrom.sizes_tobins100_50_v65_merged_color.bb
bigBedToBed hg19_nh.chrom.sizes_tobins100_50_v65_merged_color.bb clade_specific_elements_hg19.bed

# Download chain file for liftOver from hg19 to hg38
cd /home/ak2267/genomes
wget http://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz

# liftOver clade_specific_elements_hg19.bed to hg38
cd /home/ak2267/project/EnhancerAge
liftOver clade_specific_elements_hg19.bed /home/ak2267/genomes/hg19ToHg38.over.chain clade_specific_elements_hg38.bed unMapped

# Final age segmentation file is here:
/home/ak2267/project/EnhancerAge/clade_specific_elements_hg38.bed

