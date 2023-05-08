# 7/24/22
# Purpose: generate bed files for features annotated by RefSeq (use GTFs in UCSC coordinate system)
# Combine features annotated in two resources: 1) NCBI RefSeq from UCSC, and 2) UCSC RefSeq also from UCSC
# These two annotation types both come from RefSeq, but the UCSC version has been re-mapped to UCSC and may differ slightly.
# To be conservative, include both feature sets in the final feature set to be sure everything potentially a feature gets filtered.

# Download NCBI RefSeq annotations from UCSC so that chromosome names are already converted to UCSC genomes
cd /home/ak2267/genomes/RefSeq/NCBI_RefSeq
wget http://hgdownload.soe.ucsc.edu/goldenPath/rheMac10/bigZips/genes/rheMac10.ncbiRefSeq.gtf.gz
wget https://hgdownload.soe.ucsc.edu/goldenPath/calJac4/bigZips/genes/ncbiRefSeq.gtf.gz \
    ; mv ncbiRefSeq.gtf.gz calJac4.ncbiRefSeq.gtf.gz
wget https://hgdownload.soe.ucsc.edu/goldenPath/mm39/bigZips/genes/mm39.ncbiRefSeq.gtf.gz
wget https://hgdownload.soe.ucsc.edu/goldenPath/rn7/bigZips/genes/ncbiRefSeq.gtf.gz \
    ; mv ncbiRefSeq.gtf.gz rn7.ncbiRefSeq.gtf.gz
wget https://hgdownload.soe.ucsc.edu/goldenPath/susScr11/bigZips/genes/susScr11.ncbiRefSeq.gtf.gz
wget https://hgdownload.soe.ucsc.edu/goldenPath/canFam6/bigZips/genes/ncbiRefSeq.gtf.gz \
    ; mv ncbiRefSeq.gtf.gz canFam6.ncbiRefSeq.gtf.gz
wget https://hgdownload.soe.ucsc.edu/goldenPath/felCat9/bigZips/genes/felCat9.ncbiRefSeq.gtf.gz
wget https://hgdownload.soe.ucsc.edu/goldenPath/equCab3/bigZips/genes/equCab3.ncbiRefSeq.gtf.gz
wget http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/genes/hg38.ncbiRefSeq.gtf.gz

# Download UCSC RefSeq annotations from UCSC
cd /home/ak2267/genomes/RefSeq/UCSC_RefSeq
wget http://hgdownload.soe.ucsc.edu/goldenPath/rheMac10/bigZips/genes/rheMac10.refGene.gtf.gz
wget http://hgdownload.soe.ucsc.edu/goldenPath/calJac4/bigZips/genes/refGene.gtf.gz \
    ; mv refGene.gtf.gz calJac4.refGene.gtf.gz
wget http://hgdownload.soe.ucsc.edu/goldenPath/mm39/bigZips/genes/refGene.gtf.gz \
    ; mv refGene.gtf.gz mm39.refGene.gtf.gz
wget http://hgdownload.soe.ucsc.edu/goldenPath/rn7/bigZips/genes/refGene.gtf.gz \
    ; mv refGene.gtf.gz rn7.refGene.gtf.gz
wget http://hgdownload.soe.ucsc.edu/goldenPath/susScr11/bigZips/genes/susScr11.refGene.gtf.gz
wget http://hgdownload.soe.ucsc.edu/goldenPath/canFam6/bigZips/genes/refGene.gtf.gz \
    ; mv refGene.gtf.gz canFam6.refGene.gtf.gz
wget http://hgdownload.soe.ucsc.edu/goldenPath/felCat9/bigZips/genes/felCat9.refGene.gtf.gz
wget http://hgdownload.soe.ucsc.edu/goldenPath/equCab3/bigZips/genes/equCab3.refGene.gtf.gz
wget http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/genes/hg38.refGene.gtf.gz

# download these NCBI files to /home/ak2267/genomes/RefSeq/NCBI - since I downloaded them, these paths have been modified to be "refseq" instead of "refseq_old", as in human
wget https://ftp.ncbi.nih.gov/genomes/refseq_old/vertebrate_mammalian/Macaca_mulatta/latest_assembly_versions/GCF_003339765.1_Mmul_10/GCF_003339765.1_Mmul_10_genomic.gtf.gz
wget https://ftp.ncbi.nih.gov/genomes/refseq_old/vertebrate_mammalian/Callithrix_jacchus/latest_assembly_versions/GCF_009663435.1_Callithrix_jacchus_cj1700_1.1/GCF_009663435.1_Callithrix_jacchus_cj1700_1.1_genomic.gtf.gz
wget https://ftp.ncbi.nih.gov/genomes/refseq_old/vertebrate_mammalian/Mus_musculus/latest_assembly_versions/GCF_000001635.27_GRCm39/GCF_000001635.27_GRCm39_genomic.gtf.gz
wget https://ftp.ncbi.nih.gov/genomes/refseq_old/vertebrate_mammalian/Rattus_norvegicus/latest_assembly_versions/GCF_015227675.2_mRatBN7.2/GCF_015227675.2_mRatBN7.2_genomic.gtf.gz
wget https://ftp.ncbi.nih.gov/genomes/refseq_old/vertebrate_mammalian/Sus_scrofa/latest_assembly_versions/GCF_000003025.6_Sscrofa11.1/GCF_000003025.6_Sscrofa11.1_genomic.gtf.gz
wget https://ftp.ncbi.nih.gov/genomes/refseq_old/vertebrate_mammalian/Canis_lupus_familiaris/latest_assembly_versions/GCF_000002285.5_Dog10K_Boxer_Tasha/GCF_000002285.5_Dog10K_Boxer_Tasha_genomic.gtf.gz
wget https://ftp.ncbi.nih.gov/genomes/refseq_old/vertebrate_mammalian/Felis_catus/all_assembly_versions/GCF_000181335.3_Felis_catus_9.0/GCF_000181335.3_Felis_catus_9.0_genomic.gtf.gz
wget https://ftp.ncbi.nih.gov/genomes/refseq_old/vertebrate_mammalian/Equus_caballus/latest_assembly_versions/GCF_002863925.1_EquCab3.0/GCF_002863925.1_EquCab3.0_genomic.gtf.gz
wget https://ftp.ncbi.nih.gov/genomes/refseq/vertebrate_mammalian/Homo_sapiens/latest_assembly_versions/GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_genomic.gtf.gz

# STRATEGY:
# 1) loci with a 'biotype' in the GTFs from NCBI (in /home/ak2267/genomes/RefSeq/NCBI)
# --> remove exons (CDS + UTRs), 2kb upstream, plus remove introns if a pseudogene
# 2) loci without an annotation type in the GTFs from NCBI
# --> remove all of feature (i.e. including potential introns), and 2kb upstream

# run script to output blacked out regions in 6 files per genome
# ignore genes that are on patched regions since I don't know what to do with them
cd /home/ak2267/genomes/RefSeq/
python makeCategoryBedFiles_plusIntrons_chrSize.py hg38 NCBI_RefSeq/hg38.ncbiRefSeq.gtf \
       UCSC_RefSeq/hg38.refGene.gtf NCBI/GCF_000001405.40_GRCh38.p14_genomic.gtf \
       /home/ak2267/genomes/chrom.sizes/hg38.chrom.sizes
python makeCategoryBedFiles_plusIntrons_chrSize.py rheMac10 NCBI_RefSeq/rheMac10.ncbiRefSeq.gtf \
       UCSC_RefSeq/rheMac10.refGene.gtf NCBI/GCF_003339765.1_Mmul_10_genomic.gtf \
       /home/ak2267/genomes/chrom.sizes/rheMac10.chrom.sizes
python makeCategoryBedFiles_plusIntrons_chrSize.py calJac4 NCBI_RefSeq/calJac4.ncbiRefSeq.gtf \
       UCSC_RefSeq/calJac4.refGene.gtf NCBI/GCF_009663435.1_Callithrix_jacchus_cj1700_1.1_genomic.gtf \
       /home/ak2267/genomes/chrom.sizes/calJac4.chrom.sizes
python makeCategoryBedFiles_plusIntrons_chrSize.py mm39 NCBI_RefSeq/mm39.ncbiRefSeq.gtf \
       UCSC_RefSeq/mm39.refGene.gtf NCBI/GCF_000001635.27_GRCm39_genomic.gtf \
       /home/ak2267/genomes/chrom.sizes/mm39.chrom.sizes
python makeCategoryBedFiles_plusIntrons_chrSize.py rn7 NCBI_RefSeq/rn7.ncbiRefSeq.gtf \
       UCSC_RefSeq/rn7.refGene.gtf NCBI/GCF_015227675.2_mRatBN7.2_genomic.gtf \
       /home/ak2267/genomes/chrom.sizes/rn7.chrom.sizes
python makeCategoryBedFiles_plusIntrons_chrSize.py susScr11 NCBI_RefSeq/susScr11.ncbiRefSeq.gtf \
       UCSC_RefSeq/susScr11.refGene.gtf NCBI/GCF_000003025.6_Sscrofa11.1_genomic.gtf \
       /home/ak2267/genomes/chrom.sizes/susScr11.chrom.sizes
python makeCategoryBedFiles_plusIntrons_chrSize.py canFam6 NCBI_RefSeq/canFam6.ncbiRefSeq.gtf \
       UCSC_RefSeq/canFam6.refGene.gtf NCBI/GCF_000002285.5_Dog10K_Boxer_Tasha_genomic.gtf \
       /home/ak2267/genomes/chrom.sizes/canFam6.chrom.sizes
python makeCategoryBedFiles_plusIntrons_chrSize.py felCat9 NCBI_RefSeq/felCat9.ncbiRefSeq.gtf \
       UCSC_RefSeq/felCat9.refGene.gtf NCBI/GCF_000181335.3_Felis_catus_9.0_genomic.gtf \
       /home/ak2267/genomes/chrom.sizes/felCat9.chrom.sizes
python makeCategoryBedFiles_plusIntrons_chrSize.py equCab3 NCBI_RefSeq/equCab3.ncbiRefSeq.gtf \
       UCSC_RefSeq/equCab3.refGene.gtf NCBI/GCF_002863925.1_EquCab3.0_genomic.gtf \
       /home/ak2267/genomes/chrom.sizes/equCab3.chrom.sizes

# merge within each file to make final set
for species in hg38 rheMac10 calJac4 mm39 rn7 susScr11 canFam6 felCat9 equCab3
do
    echo $species
    for file in 1_proteinCoding_promoters 2_proteinCoding_exons 3_lncRNA 4_ncRNA 5_pseudogene 6_unknown
    do
	echo $file
	sort -k1,1 -k2,2n ${species}_${file}_preMerge.bed > ${species}_${file}_sort.bed
	bedtools merge -c 4 -o distinct -i ${species}_${file}_sort.bed > featureAnnotations/${species}_${file}_merged.bed
    done
done

# remove intermediate files
rm *_preMerge.bed
rm *_sort.bed

# merge 6 categories to make a single bed file for use in pipeline (to simplify bedtools commands)
# new single files are named: /home/ak2267/genomes/RefSeq/featureAnnotations/${species}_allFeatures.bed
cd /home/ak2267/genomes/RefSeq/featureAnnotations/
for species in hg38 rheMac10 calJac4 mm39 rn7 susScr11 canFam6 felCat9 equCab3
do
    echo $species
    for file in 1_proteinCoding_promoters 2_proteinCoding_exons 3_lncRNA 4_ncRNA 5_pseudogene 6_unknown
    do
	cat ${species}_${file}_merged.bed >> ${species}_allCat.bed
    done
    sort -k1,1 -k2,2n ${species}_allCat.bed > ${species}_allCat_sort.bed
    bedtools merge -c 4 -o distinct -i ${species}_allCat_sort.bed > ${species}_allFeatures.bed
done
rm *_allCat.bed
rm *_allCat_sort.bed

# convert to bigBed and upload to the browser
# added to trackDb in makeTrackDb.py

# note human is skipped (no browser tracks)
for species in rheMac10 calJac4 mm39 rn7 susScr11 canFam6 felCat9 equCab3
do
    echo $species
    # limit col 4 to 255 characters
    cut -f 1,2,3 ${species}_allFeatures.bed > temp1.bed
    cut -f 4 ${species}_allFeatures.bed | cut -c 1-255 > temp2.bed
    paste temp1.bed temp2.bed > temp3.bed
    # convert
    bedToBigBed temp3.bed /home/ak2267/genomes/chrom.sizes/${species}.chrom.sizes ${species}_allFeatures.bb
    scp ${species}_allFeatures.bb ak2267@10.5.37.220:/home/ak2267/akocher_www/CGI/${species}
done

# test if it converts back - yes they all are the same line numbers as before conversion
for species in rheMac10 calJac4 mm39 rn7 susScr11 canFam6 felCat9 equCab3
do
    bigBedToBed ${species}_allFeatures.bb ${species}_allFeatures_backToBed.bed
done


# LIFT to hg19, rheMac2, mm9
liftOver hg38_allFeatures.bed /home/ak2267/genomes/chain/hg38ToHg19.over.chain.gz hg19_allFeatures_preSort.bed unMapped
sort -k1,1 -k2,2n hg19_allFeatures_preSort.bed > hg19_allFeatures.bed
# 290381 features vs 292297 in hg38

liftOver rheMac10_allFeatures.bed /home/ak2267/genomes/chain/rheMac10ToRheMac8.over.chain.gz rheMac8_allFeatures.bed unMapped
liftOver rheMac8_allFeatures.bed /home/ak2267/genomes/chain/rheMac8ToRheMac2.over.chain.gz rheMac2_allFeatures_preSort.bed unMapped
sort -k1,1 -k2,2n rheMac2_allFeatures_preSort.bed > rheMac2_allFeatures.bed
# 209963 vs 236165 in rheMac10

liftOver mm39_allFeatures.bed /home/ak2267/genomes/chain/mm39ToMm10.over.chain.gz mm10_allFeatures.bed unMapped
liftOver mm10_allFeatures.bed /home/ak2267/genomes/chain/mm10ToMm9.over.chain.gz mm9_allFeatures_preSort.bed unMapped
sort -k1,1 -k2,2n mm9_allFeatures_preSort.bed > mm9_allFeatures.bed
# 260484 vs 261331 in mm39
