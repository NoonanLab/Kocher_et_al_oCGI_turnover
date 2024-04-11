# 3/2/2024
# Purpose: determine the proportion of gained vs lost CGIs across the tree
# for those containing a phastCons element
# This is revised manuscript Fig S21

# Pull code from 220906_CGIs_acrossEntireTree.sh where I attempted something similar
# But re-run it because since then I decided to merge oCGIs within 200 bp of each other

# GOAL: make a table with rows = locations in the human genome and columns = CGI status in all 9 species
# plus more columns with phastCons info
# then I can use this in R to count gain vs loss events in sites with / without phastCons elements

# use hg38 oCGIs and Roller genome versions

mkdir /gpfs/gibbs/pi/noonan/ak2267/Revisions/Gain_vs_Loss
cd /gpfs/gibbs/pi/noonan/ak2267/Revisions/Gain_vs_Loss

# restrict each species to intronic-intergenic CGIs only (one file for all except human & mouse, which also have FANTOM and Blacklist)
for species in rheMac10 calJac4 rn7 susScr11 canFam6 felCat9 equCab3
do
    bedtools intersect -v -a /home/ak2267/genomes/CGI/UCSC_AL/${species}_CGIsAL.bed -b /home/ak2267/genomes/RefSeq/featureAnnotations/${species}_allFeatures.bed >\
	     ${species}_CGIsAL_noFeatures.bed
done

for species in hg38 mm39
do
    bedtools intersect -v -a /home/ak2267/genomes/CGI/UCSC_AL/${species}_CGIsAL.bed \
	 -b /home/ak2267/genomes/RefSeq/featureAnnotations/${species}_allFeatures.bed \
	 -b /home/ak2267/genomes/FANTOM/FANTOM_TSS_${species}.bed \
	 -b /home/ak2267/genomes/blacklist/blacklist_${species}.bed > ${species}_CGIsAL_noFeatures.bed
done

# rename CGIs to keep track of initial species
for species in hg38 rheMac10 calJac4 mm39 rn7 susScr11 canFam6 felCat9 equCab3
do
    echo $species
    awk '{ print $1"\t"$2"\t"$3"\t'$species'_"NR }' ${species}_CGIsAL_noFeatures.bed > ${species}_CGIsAL_named.bed
done
rm *_noFeatures.bed

# lift to human (hg38)
for species in rheMac10 calJac4 mm39 rn7 susScr11 canFam6 felCat9 equCab3
do
    echo ${species}
    liftOver -minMatch=0.3 ${species}_CGIsAL_named.bed /home/ak2267/genomes/chain/${species}ToHg38.over.chain.gz ${species}_CGIsAL_named_LOhg38.bed unMapped
done

# lift back to initial species, intersect with original sites, and filter to only those that do (in hg38 coordinates)
# makes files ${species}_CGIsThatMapBack_inHg38.bed
module load Python/3.8.6-GCCcore-10.2.0

for species in rheMac10 calJac4 mm39 rn7 susScr11 canFam6 felCat9 equCab3
do 
    echo $species
    liftOver -minMatch=0.3 ${species}_CGIsAL_named_LOhg38.bed /home/ak2267/genomes/chain/hg38To${species^}.over.chain.gz ${species}_CGIsAL_named_LOhg38_liftBack.bed unMapped
    bedtools intersect -wao -a ${species}_CGIsAL_named_LOhg38_liftBack.bed -b ${species}_CGIsAL_named.bed > ${species}_checkMappingBack.txt
    python /gpfs/gibbs/pi/noonan/ak2267/speciesPairs/consensusCGIs/Scripts/restrictToSitesThatMapBack.py ${species}_checkMappingBack.txt > ${species}_CGIsThatMapBack_inOwnCoord.bed
    python /gpfs/gibbs/pi/noonan/ak2267/speciesPairs/consensusCGIs/Scripts/restrictToLO.py ${species}_CGIsAL_named_LOhg38.bed ${species}_CGIsThatMapBack_inOwnCoord.bed > ${species}_CGIsThatMapBack_inHg38.bed
done

# merge in human, turning names into comma-separated list to keep track of species
for species in rheMac10 calJac4 mm39 rn7 susScr11 canFam6 felCat9 equCab3
do
    cat ${species}_CGIsThatMapBack_inHg38.bed >> cat_withoutHuman.bed
done
cat cat_withoutHuman.bed hg38_CGIsAL_named.bed > cat.bed # add human
sort -k1,1 -k2,2n cat.bed > sort.bed
bedtools merge -c 4 -o distinct -i sort.bed > mergedCGIs_hg38.bed ##### THIS IS ALL CGIs IN hg38 COORD

# lift out to each original species to determine which have the sequence present
for species in rheMac10 calJac4 mm39 rn7 susScr11 canFam6 felCat9 equCab3
do
    echo $species
    liftOver -minMatch=0.3 mergedCGIs_hg38.bed /home/ak2267/genomes/chain/hg38To${species^}.over.chain.gz mergedCGIs_hg38_LO${species}.bed unMapped
done

# lift back to human and filter based on whether it intersects the original site - submit as jobs because these files are large
for species in rheMac10 calJac4 mm39 rn7 susScr11 canFam6 felCat9 equCab3
do
    echo 'cd /gpfs/gibbs/pi/noonan/ak2267/Revisions/Gain_vs_Loss ; source /home/ak2267/.bashrc ; source /home/ak2267/.bash_profile ; liftOver -minMatch=0.3 mergedCGIs_hg38_LO'${species}'.bed /home/ak2267/genomes/chain/'${species}'ToHg38.over.chain.gz mergedCGIs_hg38_LO'${species}'_liftBackToHg38.bed unMapped'
done >> 240302_liftBackToHuman_jobFile.txt
dsq --job-file 240302_liftBackToHuman_jobFile.txt --mem-per-cpu 10G -c 1 --mail-type FAIL,END
sbatch sbatch dsq-240302_liftBackToHuman_jobFile-2024-03-02.sh # 22669673

for species in rheMac10 calJac4 mm39 rn7 susScr11 canFam6 felCat9 equCab3
do
    echo 'cd /gpfs/gibbs/pi/noonan/ak2267/Revisions/Gain_vs_Loss ; source /home/ak2267/.bashrc ; source /home/ak2267/.bash_profile ; bedtools intersect -wao -a mergedCGIs_hg38_LO'${species}'_liftBackToHg38.bed -b mergedCGIs_hg38.bed > '${species}'_checkMappingBack_human.txt ; python /gpfs/gibbs/pi/noonan/ak2267/speciesPairs/consensusCGIs/Scripts/restrictToSitesThatMapBack.py '${species}'_checkMappingBack_human.txt > '${species}'_CGIsThatMapBack_inHumanCoord.bed ; python /gpfs/gibbs/pi/noonan/ak2267/speciesPairs/consensusCGIs/Scripts/restrictToLO.py mergedCGIs_hg38_LO'${species}'.bed '${species}'_CGIsThatMapBack_inHumanCoord.bed > '${species}'_CGIsThatMapBackToHuman_inOwnCoord.bed'
done >> 240302_restrictInHuman_jobFile.txt
dsq --job-file 240302_restrictInHuman_jobFile.txt --mem-per-cpu 5G -c 1 --mail-type FAIL,END
sbatch dsq-240302_restrictInHuman_jobFile-2024-03-02.sh # 22671589


### filter based on not having a feature in any of the species
bedtools intersect -v -a mm39_CGIsThatMapBackToHuman_inOwnCoord.bed \
	 -b /home/ak2267/genomes/RefSeq/featureAnnotations/mm39_allFeatures.bed \
	 -b /home/ak2267/genomes/FANTOM/FANTOM_TSS_mm39.bed \
	 -b /home/ak2267/genomes/blacklist/blacklist_mm39.bed > mm39_mergedCGIs_noFeatures_inmm39.bed
bedtools intersect -v -a mergedCGIs_hg38.bed \
	 -b /home/ak2267/genomes/RefSeq/featureAnnotations/hg38_allFeatures.bed \
	 -b /home/ak2267/genomes/FANTOM/FANTOM_TSS_hg38.bed \
	 -b /home/ak2267/genomes/blacklist/blacklist_hg38.bed > hg38_mergedCGIs_noFeatures_inhg38.bed
for species in rheMac10 calJac4 rn7 susScr11 canFam6 felCat9 equCab3
do
    bedtools intersect -v -a ${species}_CGIsThatMapBackToHuman_inOwnCoord.bed \
	 -b /home/ak2267/genomes/RefSeq/featureAnnotations/${species}_allFeatures.bed > ${species}_mergedCGIs_noFeatures_in${species}.bed
done

# WHAT WE HAVE CREATED:
# CGIs in human coordinates: mergedCGIs_hg38.bed
# CGIs that lift to each species: ${species}_CGIsThatMapBackToHuman_inOwnCoord.bed
# CGIs that lift to each species and don't overlap a feature in those coordinates: ${species}_mergedCGIs_noFeatures_in${species}.bed


# get phastCons info for each CGI
bedtools intersect -wa -u -a mergedCGIs_hg38.bed -b /home/ak2267/genomes/phastCons/phastConsElements100way_hg38.bed > mergedCGIs_hg38_overlapPhastCons.bed


# make test file
grep 'chr12' mergedCGIs_hg38.bed > mergedCGIs_chr12_hg38.bed

# use a python script to summarize across entire tree. makes table with:
# rows = CGIs (names from col4)
# columns = species, entry in each box = 0 vs 1 vs 2 (no sequence, sequence but no CGI, sequence with CGI)
# exclude CGIs if they overlap a feature in any species

# echo 'cd /gpfs/gibbs/pi/noonan/ak2267/Revisions/Gain_vs_Loss ; source /home/ak2267/.bashrc ; source /home/ak2267/.bash_profile ; python ak20240302_summarizeCGIsAcrossEntireTree.py rheMac10,calJac4,mm39,rn7,susScr11,canFam6,felCat9,equCab3' >> 240302_summaryTable_jobFile.txt
# dsq --job-file 240302_summaryTable_jobFile.txt --mem-per-cpu 10G -c 1 --mail-type FAIL,END --time=03:00:00
# sbatch dsq-240302_summaryTable_jobFile-2024-03-03.sh # 

# I ran this on an interactive node in <10 min
python ak20240302_summarizeCGIsAcrossEntireTree.py rheMac10,calJac4,mm39,rn7,susScr11,canFam6,felCat9,equCab3

# download output file for plotting in R
cd /Users/acadiak/Desktop/Yale/\!Writing/Kocher_CpG/GB_Revisions/Revised_code/Step12_revisions/Gain_vs_Loss
scp ak2267@mccleary.ycrc.yale.edu:/gpfs/gibbs/pi/noonan/ak2267/Revisions/Gain_vs_Loss/ak20240302_CGI_summary_acrossTree.txt .

# go into R for further plotting
# ak20240302_gain_vs_loss.R