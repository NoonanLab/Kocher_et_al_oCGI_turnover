# 7/25/22
# Purpose: download RepeatMasker files (all genomes) for use in downstream CGI analysis

### RepeatMasker

# Download rmsk.txt files
mkdir /home/ak2267/genomes/rmsk
cd /home/ak2267/genomes/rmsk

for species in rheMac10 calJac4 mm39 rn7 susScr11 canFam6 felCat9 equCab3 hg38
do
    wget https://hgdownload.soe.ucsc.edu/goldenPath/${species}/database/rmsk.txt.gz
    mv rmsk.txt.gz rmsk_${species}.txt.gz
done

# Convert to bed files
for species in rheMac10 calJac4 mm39 rn7 susScr11 canFam6 felCat9 equCab3 hg38
do
    echo ${species}
    zcat rmsk_${species}.txt.gz | awk '{ print $6"\t"$7"\t"$8"\t"$12";"$13";"$11 }' > rmsk_${species}.bed
done

# path for these files: /home/ak2267/genomes/rmsk/rmsk_${species}.bed

# download for Noonan genomes on 7/29/22
for species in hg19 rheMac2 mm10
do
    wget https://hgdownload.soe.ucsc.edu/goldenPath/${species}/database/rmsk.txt.gz
    mv rmsk.txt.gz rmsk_${species}.txt.gz
done

# Convert to bed files
for species in hg19 rheMac2 mm10
do
    echo ${species}
    zcat rmsk_${species}.txt.gz | awk '{ print $6"\t"$7"\t"$8"\t"$12";"$13";"$11 }' > rmsk_${species}.bed
done

# 8/4/22
# merge entries so that there are no overlapping regions - necessary for calculating % of CGI covered by repeats later on
cd /home/ak2267/genomes/rmsk
for species in rheMac10 calJac4 mm39 rn7 susScr11 canFam6 felCat9 equCab3 hg38 hg19 rheMac2 mm10
do
    echo ${species}
    bedtools merge -c 4 -o collapse -i rmsk_${species}.bed > rmsk_${species}_merged.bed
done

# lift mm10 to mm9 (only mm10 was available from UCSC)
liftOver rmsk_mm10.bed /home/ak2267/genomes/chain/mm10ToMm9.over.chain.gz rmsk_mm9.bed unMapped
liftOver rmsk_mm10_merged.bed /home/ak2267/genomes/chain/mm10ToMm9.over.chain.gz rmsk_mm9_merged.bed unMapped


