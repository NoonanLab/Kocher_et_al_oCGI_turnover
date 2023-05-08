# 7/25/22
# Purpose: download ENCODE blacklist (human and mouse) and format for hg38/hg19 and mm39/mm9

### ENCODE blacklist
# See here for description
# https://sites.google.com/site/anshulkundaje/projects/blacklists

mkdir /home/ak2267/genomes/blacklist
cd /home/ak2267/genomes/blacklist

# Download human (hg38)
# https://www.encodeproject.org/files/ENCFF356LFX/
wget https://www.encodeproject.org/files/ENCFF356LFX/@@download/ENCFF356LFX.bed.gz
mv ENCFF356LFX.bed.gz blacklist_hg38.bed.gz

# Download mouse (mm10)
# https://github.com/Boyle-Lab/Blacklist/tree/master/lists
wget https://github.com/Boyle-Lab/Blacklist/raw/master/lists/mm10-blacklist.v2.bed.gz
mv mm10-blacklist.v2.bed.gz blacklist_mm10.bed.gz

# Lift to all other required species (hg38 --> hg19) and (mm10 --> mm39, mm9)
gunzip *.gz
liftOver blacklist_hg38.bed /home/ak2267/genomes/chain/hg38ToHg19.over.chain.gz blacklist_hg19.bed unMapped
# 910 --> 836 regions
cut -f 1,2,3 blacklist_mm10.bed > blacklist_mm10_3col.bed
liftOver blacklist_mm10_3col.bed /home/ak2267/genomes/chain/mm10ToMm39.over.chain.gz blacklist_mm39.bed unMapped
# 3435 --> 3360 regions
liftOver blacklist_mm10_3col.bed /home/ak2267/genomes/chain/mm10ToMm9.over.chain.gz blacklist_mm9.bed unMapped
# 3435 --> 3336 regions
