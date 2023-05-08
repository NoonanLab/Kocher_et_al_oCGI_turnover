# 6/16/22
# Purpose: download FANTOM data and process for use in:
# excluding unannotated promoters from analysis

# work here
cd /home/ak2267/genomes/FANTOM

# FILES TO USE:
TSS_human_hg38.bed
TSS_mouse_mm39.bed

# Promoter data (human and mouse)
# https://fantom.gsc.riken.jp/5/datafiles/phase1.3/extra/TSS_classifier/

# this is in hg19
wget https://fantom.gsc.riken.jp/5/datafiles/phase1.3/extra/TSS_classifier/TSS_human.bed.gz
gunzip *.gz

# lift to hg38
tail -1048124 TSS_human.bed > TSS_human_noHeader.bed
liftOver TSS_human_noHeader.bed ~/genomes/chain/hg19ToHg38.over.chain.gz TSS_human_hg38.bed unMapped

# this is in mm9
wget https://fantom.gsc.riken.jp/5/datafiles/phase1.3/extra/TSS_classifier/TSS_mouse.bed.gz
gunzip *.gz

# lift to mm39
liftOver TSS_mouse.bed ~/genomes/chain/mm9ToMm10.over.chain.gz TSS_mouse_mm10.bed unMapped
liftOver TSS_mouse_mm10.bed ~/genomes/chain/mm10ToMm39.over.chain.gz TSS_mouse_mm39.bed unMapped


###### 7/25/22
# convert TSS files to hg38, hg19 and mm39, mm9 with standard name format: FANTOM_TSS_${species}.bed
sort -k1,1 -k2,2n TSS_human_noHeader.bed > FANTOM_TSS_hg19.bed # 1048124 sites
sort -k1,1 -k2,2n TSS_human_hg38.bed > FANTOM_TSS_hg38.bed # 1047890 sites
sort -k1,1 -k2,2n TSS_mouse.bed > FANTOM_TSS_mm9.bed # 652852 sites
sort -k1,1 -k2,2n TSS_mouse_mm39.bed > FANTOM_TSS_mm39.bed # 652803 sites


