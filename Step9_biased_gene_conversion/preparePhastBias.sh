# 9/14/22
# Purpose: prepare phastBias tracts for further analysis

# possible note for future: phastBias only picks up 25-50% of bases affected by gBGC
# see https://genome-asia.ucsc.edu/cgi-bin/hgTrackUi?hgsid=754326795_M3Zx1R8sS3SlqpVX6ZU5SlNJN6A3&db=hg19&c=chrX&g=phastBias

# work here
cd /gpfs/gibbs/pi/noonan/ak2267/speciesPairs/SequenceAnalysis/gBGC_analysis/

# phastBias files are here
# /home/su57/scratch60/200120_120mammals/220912_phastBias/

# turn gff files into bed files
mkdir phastBias_bed
rm phastBias_bed/*
for species in hg38 rheMac8 calJac3 mm10 rn6 susScr11 canFam3 felCat8 HLequCab3
do
    for num in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X
    do
	awk '{ print "chr'${num}'\t"$4"\t"$5"\t"$6 }' /home/su57/scratch60/200120_120mammals/220912_phastBias/220913_chr${num}_${species}_phastBias.gff > phastBias_bed/${species}_chr${num}_phastBias.bed
	cat phastBias_bed/${species}_chr${num}_phastBias.bed >> phastBias_bed/${species}_allChr_phastBias.bed
    done
done
rm phastBias_bed/*_chr*_phastBias.bed

# merge within 1kb
for species in hg38 rheMac8 calJac3 mm10 rn6 susScr11 canFam3 felCat8 HLequCab3
do
    bedtools merge -c 4 -o collapse -d 1000 -i phastBias_bed/${species}_allChr_phastBias.bed > phastBias_bed/${species}_allChr_phastBias_merged.bed
done

# lift to hg19 for use with Noonan data - note that 49/5490 fail to lift for human and 138/29369 for rhesus
cd /gpfs/gibbs/pi/noonan/ak2267/speciesPairs/SequenceAnalysis/gBGC_analysis/phastBias_bed/
liftOver -minMatch=0.8 hg38_allChr_phastBias_merged.bed /home/ak2267/genomes/chain/hg38ToHg19.over.chain.gz hg38inhg19_allChr_phastBias_merged.bed unMapped
liftOver -minMatch=0.8 rheMac8_allChr_phastBias_merged.bed /home/ak2267/genomes/chain/hg38ToHg19.over.chain.gz rheMac8inhg19_allChr_phastBias_merged.bed unMapped

# FINAL FILES (NOTE THAT ALL ARE IN hg38 COORDINATES)
/gpfs/gibbs/pi/noonan/ak2267/speciesPairs/SequenceAnalysis/gBGC_analysis/phastBias_bed/${species}_allChr_phastBias_merged.bed


