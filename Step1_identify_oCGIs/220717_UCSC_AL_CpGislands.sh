# 7/17/22
# Purpose: use Andy Law's (AL) script from UCSC to permissively call CpG islands
# using criteria from Gardiner-Garden & Frommer (1987)
# Instructions from google group question are here
# https://groups.google.com/a/soe.ucsc.edu/g/genome/c/HoERMqdfe-g
# Example code from running on galGal3 is here
# https://github.com/ucscGenomeBrowser/kent/blob/master/src/hg/makeDb/doc/galGal3.txt#L649

cd /gpfs/gibbs/pi/noonan/ak2267/genomes/
# get masked fasta files for all genomes (start with Roller genome versions plus hg38)
for genome in hg38 rheMac10 calJac4 mm39 rn7 susScr11 canFam6 felCat9 equCab3
do
    wget https://hgdownload.soe.ucsc.edu/goldenPath/${genome}/bigZips/${genome}.fa.masked.gz
done
gunzip *.fa.masked.gz

# download program to format fasta files
cd /home/ak2267/bin
wget https://hgwdev-gperez2.gi.ucsc.edu/~gperez2/mlq/mlq_29758/preProcGgfAndy.linux_x86_64
chmod +x preProcGgfAndy.linux_x86_64

# download AL script 
cd /home/ak2267/genomes/CGI/UCSC_AL
wget https://hgwdev-gperez2.gi.ucsc.edu/~gperez2/mlq/mlq_29758/ggf-andy-cpg-island.pl
chmod +x ggf-andy-cpg-island.pl

# split masked fasta into individual chromosome files
cd /home/ak2267/genomes/CGI/UCSC_AL
for species in rheMac10 calJac4 mm39 rn7 susScr11 canFam6 felCat9 equCab3 hg38
do
    echo $species
    mkdir ${species}
    faSplit byname /gpfs/gibbs/pi/noonan/ak2267/genomes/${species}.fa.masked ${species}/
done

# run UCSC prep file (preProcGgfAndy) and perl script to call CpG islands (ggf-andy-cpg-island.pl)
# had to remove line breaks from perl one-liner and then this worked great
cd /home/ak2267/genomes/CGI/UCSC_AL
module load Perl/5.32.0-GCCcore-10.2.0
for species in rheMac10 calJac4 mm39 rn7 susScr11 canFam6 felCat9 equCab3 hg38
do
    echo $species
    for file in $(ls ${species}/chr*.fa)
    do
	echo $file
	chr=$(basename $file .fa)
	/home/ak2267/bin/preProcGgfAndy.linux_x86_64 $file | ./ggf-andy-cpg-island.pl \
	    | perl -wpe 'chomp; ($s,$e,$cpg,$n,$c,$g,$oE) = split("\t"); $s--; $gc = $c + $g;  $pCpG = (100.0 * 2 * $cpg / $n); $pGc = (100.0 * $gc / $n); $_ = "'$chr'\t$s\t$e\tCpG: $cpg\t$n\t$cpg\t$gc\t" . "$pCpG\t$pGc\t$oE\n";' >> ${species}_CGIsAL_unsorted.bed
    done;done

# sort bed files
for species in rheMac10 calJac4 mm39 rn7 susScr11 canFam6 felCat9 equCab3 hg38
do
    echo $species
    sort -k1,1 -k2,2n ${species}_CGIsAL_unsorted.bed > ${species}_CGIsAL.bed
done

# merge within 200 bp
for species in rheMac10 calJac4 mm39 rn7 susScr11 canFam6 felCat9 equCab3 hg38
do
    echo $species
    bedtools merge -d 200 -i ${species}_CGIsAL.bed > ${species}_CGIsAL_merged200.bed
done

# show difference after merging
for species in rheMac10 calJac4 mm39 rn7 susScr11 canFam6 felCat9 equCab3 hg38
do
    noMerge=$(cat ${species}_CGIsAL.bed | wc -l)
    yesMerge=$(cat ${species}_CGIsAL_merged200.bed | wc -l)
    echo $species $noMerge $yesMerge
done

# rheMac10 89562 77529
# calJac4 82827 71899
# mm39 65548 58902
# rn7 91858 83399
# susScr11 161000 134252
# canFam6 127594 104994
# felCat9 218352 180238
# equCab3 153498 133232
# hg38 88722 75980

# rename merged set as regular names
for species in rheMac10 calJac4 mm39 rn7 susScr11 canFam6 felCat9 equCab3 hg38
do
    mv ${species}_CGIsAL_merged200.bed ${species}_CGIsAL.bed
done

# convert bed to bigBed for viewing on the genome browser
for species in rheMac10 calJac4 mm39 rn7 susScr11 canFam6 felCat9 equCab3 hg38
do
    echo $species
    cut -f 1,2,3,6 ${species}_CGIsAL.bed > temp_${species}_CGIsAL.bed
    bedToBigBed temp_${species}_CGIsAL.bed /home/ak2267/genomes/chrom.sizes/${species}.chrom.sizes ${species}_CGIsAL.bb
done

# transfer to lab server for viewing on the browser - NOTE I didn't do hg38
for species in rheMac10 calJac4 mm39 rn7 susScr11 canFam6 felCat9 equCab3
do
    echo $species
    scp ${species}_CGIsAL.bb ak2267@10.5.37.220:/home/ak2267/akocher_www/CGI/${species}
done

######################################################################
## call CGIs again in hg19, rheMac2, and mm9 for use with Reilly/Cotney data which has already been processed in these genomes

# get masked fasta for hg19, rheMac2, and mm9
# goals:
# 1) single species.fa.masked file in /gpfs/gibbs/pi/noonan/ak2267/genomes/
# 2) directory species/ with individual chromosome fasta files in /home/ak2267/genomes/CGI/UCSC_AL/
cd /gpfs/gibbs/pi/noonan/ak2267/genomes/

# hg19
cd /gpfs/gibbs/pi/noonan/ak2267/genomes/
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.masked.gz
gunzip hg19.fa.masked.gz
mkdir /home/ak2267/genomes/CGI/UCSC_AL/hg19
cd /home/ak2267/genomes/CGI/UCSC_AL
mkdir /home/ak2267/genomes/CGI/UCSC_AL/hg19
faSplit byname /gpfs/gibbs/pi/noonan/ak2267/genomes/hg19.fa.masked hg19/

# rheMac2
cd /gpfs/gibbs/pi/noonan/ak2267/genomes/
mkdir /home/ak2267/genomes/CGI/UCSC_AL/rheMac2
wget https://hgdownload.soe.ucsc.edu/goldenPath/rheMac2/bigZips/chromFaMasked.tar.gz
tar -xvzf chromFaMasked.tar.gz
for file in $(ls hardMask/*)
do
    fileStem=`basename $file .fa.masked`
    mv hardMask/${fileStem}.fa.masked /home/ak2267/genomes/CGI/UCSC_AL/rheMac2/${fileStem}.fa
done
rm chromFaMasked.tar.gz
for i in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 Ur X
do
    cat /home/ak2267/genomes/CGI/UCSC_AL/rheMac2/chr${i}.fa >> rheMac2.fa.masked
done

# mm9
cd /gpfs/gibbs/pi/noonan/ak2267/genomes/
mkdir /home/ak2267/genomes/CGI/UCSC_AL/mm9
wget https://hgdownload.soe.ucsc.edu/goldenPath/mm9/bigZips/chromFaMasked.tar.gz
tar -xvzf chromFaMasked.tar.gz
for file in $(ls chr*.fa.masked)
do
    fileStem=`basename $file .fa.masked`
    mv ${fileStem}.fa.masked /home/ak2267/genomes/CGI/UCSC_AL/mm9/${fileStem}.fa
done
rm chromFaMasked.tar.gz
for file in $(ls /home/ak2267/genomes/CGI/UCSC_AL/mm9/*.fa)
do
    cat ${file} >> mm9.fa.masked
done


# run UCSC prep file (preProcGgfAndy) and perl script to call CpG islands (ggf-andy-cpg-island.pl)
# had to remove line breaks from perl one-liner and then this worked great
cd /home/ak2267/genomes/CGI/UCSC_AL
module load Perl/5.32.0-GCCcore-10.2.0
for species in hg19 rheMac2 mm9
do
    echo $species
    for file in $(ls ${species}/chr*.fa)
    do
	echo $file
	chr=$(basename $file .fa)
	/home/ak2267/bin/preProcGgfAndy.linux_x86_64 $file | ./ggf-andy-cpg-island.pl \
	    | perl -wpe 'chomp; ($s,$e,$cpg,$n,$c,$g,$oE) = split("\t"); $s--; $gc = $c + $g;  $pCpG = (100.0 * 2 * $cpg / $n); $pGc = (100.0 * $gc / $n); $_ = "'$chr'\t$s\t$e\tCpG: $cpg\t$n\t$cpg\t$gc\t" . "$pCpG\t$pGc\t$oE\n";' >> ${species}_CGIsAL_unsorted.bed
    done;done

# sort bed files
for species in hg19 rheMac2 mm9
do
    echo $species
    sort -k1,1 -k2,2n ${species}_CGIsAL_unsorted.bed > ${species}_CGIsAL.bed
done

# merge within 200 bp
for species in hg19 rheMac2 mm9
do
    echo $species
    bedtools merge -d 200 -i ${species}_CGIsAL.bed > ${species}_CGIsAL_merged200.bed
done

# show difference after merging
for species in hg19 rheMac2 mm9
do
    noMerge=$(cat ${species}_CGIsAL.bed | wc -l)
    yesMerge=$(cat ${species}_CGIsAL_merged200.bed | wc -l)
    echo $species $noMerge $yesMerge
done

# hg19 81794 69956
# rheMac2 77834 68161
# mm9 67102 60183

# rename merged set as regular names
for species in hg19 rheMac2 mm9
do
    mv ${species}_CGIsAL_merged200.bed ${species}_CGIsAL.bed
done


########
# compare numbers from calling de novo in hg19/rheMac2/mm9 to simply lifting from hg38/rheMac10/mm9
# note that this was done with the UNmerged set of CGIs

cd /home/ak2267/genomes/CGI/UCSC_AL

cut -f 1,2,3,6 hg38_CGIsAL.bed > temp_hg38.bed
liftOver temp_hg38.bed /home/ak2267/genomes/chain/hg38ToHg19.over.chain.gz hg38_LOhg19_CGIsAL_preSort.bed unMapped
sort -k1,1 -k2,2n hg38_LOhg19_CGIsAL_preSort.bed > hg38_LOhg19_CGIsAL.bed
wc -l hg38_LOhg19_CGIsAL.bed # 87493 lifted vs 81794 de novo in hg19

cut -f 1,2,3,6 rheMac10_CGIsAL.bed > temp_rheMac10.bed
liftOver temp_rheMac10.bed /home/ak2267/genomes/chain/rheMac10ToRheMac8.over.chain.gz temp_rheMac10_LOrheMac8.bed unMapped
liftOver temp_rheMac10_LOrheMac8.bed /home/ak2267/genomes/chain/rheMac8ToRheMac2.over.chain.gz rheMac10_LOrheMac2_preSort.bed unMapped
sort -k1,1 -k2,2n rheMac10_LOrheMac2_preSort.bed > rheMac10_LOrheMac2_CGIsAL.bed
wc -l rheMac10_LOrheMac2_CGIsAL.bed # 71285 lifted vs 77834 de novo in rheMac2

cut -f 1,2,3,6 mm39_CGIsAL.bed > temp_mm39.bed
liftOver temp_mm39.bed /home/ak2267/genomes/chain/mm39ToMm10.over.chain.gz temp_mm39_LOmm10.bed unMapped
liftOver temp_mm39_LOmm10.bed /home/ak2267/genomes/chain/mm10ToMm9.over.chain.gz mm39_LOmm9_preSort.bed unMapped
sort -k1,1 -k2,2n mm39_LOmm9_preSort.bed > mm39_LOmm9_CGIsAL.bed
wc -l mm39_LOmm9_CGIsAL.bed # 65368 lifted vs 67102 de novo in rheMac2

# RESULT: more from lifting in human, but more de novo in rhesus and mouse
# DECISION: use de novo set!
/home/ak2267/genomes/CGI/UCSC_AL/hg19_CGIsAL.bed
/home/ak2267/genomes/CGI/UCSC_AL/rheMac2_CGIsAL.bed
/home/ak2267/genomes/CGI/UCSC_AL/mm9_CGIsAL.bed

# 8/10/22
# get unmasked fasta files for use downstream

cd /gpfs/gibbs/pi/noonan/ak2267/genomes
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz
wget https://hgdownload.soe.ucsc.edu/goldenPath/rheMac2/bigZips/rheMac2.fa.gz
wget https://hgdownload.soe.ucsc.edu/goldenPath/mm9/bigZips/mm9.fa.gz
gunzip *.gz
