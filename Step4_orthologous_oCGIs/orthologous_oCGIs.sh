# 7/30/22
# Purpose: identify orthologous CGIs between all species pairs (called "consensusCGIs" in this pipeline)
# This means the 8 Roller species (lifting via human, hg38): rheMac10 calJac4 mm39 rn7 susScr11 canFam6 felCat9 equCab3
# And the 3 species for Noonan data (lifting via human, hg19): hg19 rheMac2 mm9

# NOTE that this pipeline does NOT consider some things that need to be incorporated downstream:
# removing CGIs that overlap a promoter peak (want to do separately for each mark)
# annotating with rmsk overlap in each species

#### work here
mkdir /gpfs/gibbs/pi/noonan/ak2267/speciesPairs
mkdir /gpfs/gibbs/pi/noonan/ak2267/speciesPairs/consensusCGIs
cd /gpfs/gibbs/pi/noonan/ak2267/speciesPairs/consensusCGIs

#### files to incorporate
# CGIs (AL version)
/home/ak2267/genomes/CGI/UCSC_AL/${species}_CGIsAL.bed
# RefSeq featureAnnotations
/home/ak2267/genomes/RefSeq/featureAnnotations/${species}_allFeatures.bed
# FANTOM TSS annotation
/home/ak2267/genomes/FANTOM/FANTOM_TSS_${species}.bed
# ENCODE blacklist
/home/ak2267/genomes/blacklist/blacklist_${species}.bed

#### run prepFiles_consensusCGIs.py to make job file
# will have 28 jobs for Roller genomes, plus 3 jobs for Noonan genomes
cd /gpfs/gibbs/pi/noonan/ak2267/speciesPairs/consensusCGIs

# with merge of CGIs within 200bp
python prepFiles_consensusCGIs.py > 220803_consensusCGIs_jobFile_200.txt
dsq --job-file 220803_consensusCGIs_jobFile_200.txt --mem-per-cpu 5G -c 1 --mail-type FAIL,END
sbatch dsq-220803_consensusCGIs_jobFile_200-2023-02-11.sh # 21095526 - run with NEW version of restrict script

# this makes files for downstream use with the names:
/gpfs/gibbs/pi/noonan/ak2267/speciesPairs/consensusCGIs/{speciesA}_{speciesB}/{speciesA}_{speciesB}_speciesA_reconciledCGI.bed
/gpfs/gibbs/pi/noonan/ak2267/speciesPairs/consensusCGIs/{speciesA}_{speciesB}/{speciesA}_{speciesB}_speciesB_reconciledCGI.bed
/gpfs/gibbs/pi/noonan/ak2267/speciesPairs/consensusCGIs/{speciesA}_{speciesB}/{speciesA}_{speciesB}_{hg38/hg19}_reconciledCGI.bed # --> for use with phastCons elements

