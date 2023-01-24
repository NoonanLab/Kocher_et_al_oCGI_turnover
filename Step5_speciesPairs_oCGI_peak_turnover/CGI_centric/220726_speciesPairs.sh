# 7/26/22
# Purpose: make summary tables for each species pair and the accompanying ChIP data in each tissue
# Result will be one summary table for each species pair x tissue x mark (and x time point for Noonan data)

# Orthologous (aka "reconciled") CGIs are here:
/gpfs/gibbs/pi/noonan/ak2267/speciesPairs/consensusCGIs/{speciesA}_{speciesB}/{speciesA}_{speciesB}_speciesA_reconciledCGI.bed
/gpfs/gibbs/pi/noonan/ak2267/speciesPairs/consensusCGIs/{speciesA}_{speciesB}/{speciesA}_{speciesB}_speciesB_reconciledCGI.bed

################ CGI-centric comparisons
# Work here:
mkdir /gpfs/gibbs/pi/noonan/ak2267/speciesPairs/CGIcentric
cd /gpfs/gibbs/pi/noonan/ak2267/speciesPairs/CGIcentric

#### run prepFiles_speciesPairs_CGIcentric.py
# Will have
# Roller: 28 pairs x 4 tissues x 3 marks = 336 jobs
# Noonan - Reilly: 1 pair x 4 tissues x 2 marks + 2 pairs x 3 tissues x 2 marks = 20 jobs
# Noonan - Cotney: 1 pair x 4 timePoints x 1 mark + 2 pairs x 2 timePoints x 1 mark = 8 jobs
# CTCF: 6 pairs x 1 tissue x 1 TF = 6 jobs
# Liver TFs: 6 pairs x 1 tissue x 4 TFs = 24 jobs

# Roller ChIP peaks and signal:
/gpfs/gibbs/pi/noonan/ak2267/Roller/peaks/intersection/${species}_${tissue}_${mark}_intersection.bed
/gpfs/gibbs/pi/noonan/ak2267/Roller/peaks/repsToUseRoller.txt
/gpfs/gibbs/pi/noonan/ak2267/Roller/bam/${species}_${tissue}_${mark}_${replicate}_bowtie2_filtered.bam
# CTCF ChIP peaks and signal:
/gpfs/gibbs/pi/noonan/ak2267/LiverTF/peaks/CTCF/intersection/${species}_liver_CTCF_intersection.bed
/gpfs/gibbs/pi/noonan/ak2267/LiverTF/bam/CTCF/${species}_liver_CTCF_${replicate}_filtered.bam
# Liver TF ChIP peaks and signal:
/gpfs/gibbs/pi/noonan/ak2267/LiverTF/peaks/intersection/${species}_liver_${TF}_intersection.bed
/gpfs/gibbs/pi/noonan/ak2267/LiverTF/peaks/repsToUseLiverTF.txt
/gpfs/gibbs/pi/noonan/ak2267/LiverTF/bam/${species}_${TF}_${replicate}_filtered.bam
# Noonan BRAIN ChIP peaks and signal:
/home/ak2267/project/EnhancerClasses/hg19/ac/merge_${timePoint}_overlap_named.bed # timePoint = CS16/CS23/F2F/F2O
/home/ak2267/project/EnhancerClasses/hg19/me2/merge_${timePoint}_me2_overlap_named.bed # timePoint = CS16/CS23/F2F/F2O
/home/ak2267/project/EnhancerClasses/rheMac2/ac/merge_${timePoint}_overlap_rh_named.bed # timePoint = e55/e79F/e79O
/home/ak2267/project/EnhancerClasses/rheMac2/me2/merge_${timePoint}_overlap_me2_rh_named.bed # timePoint = e55/e79F/e79O
/home/ak2267/project/EnhancerClasses/mm9/ac/merge_${timePoint}_overlap_mm_named.bed # timePoint = e11/e14/17F/17O
/home/ak2267/project/EnhancerClasses/mm9/me2/merge_${timePoint}_overlap_me2_mm_named.bed # timePoint = e11/e14/17F/17O
/gpfs/gibbs/pi/noonan/ak2267/NoonanBrain/hg19/ac/${timePoint}_ac_rep{1/2}.bw # timePoint = CS16/CS23/F2F/F2O
/gpfs/gibbs/pi/noonan/ak2267/NoonanBrain/hg19/me2/${timePoint}_me2_rep{1/2}.bw # timePoint = CS16/CS23/F2F/F2O
/gpfs/gibbs/pi/noonan/ak2267/NoonanBrain/rheMac2/ac/${timePoint}_ac_rep{1/2}.bw # timePoint = e55/e79F/e79O - only rep1 for e55
/gpfs/gibbs/pi/noonan/ak2267/NoonanBrain/rheMac2/me2/${timePoint}_me2_rep{1/2}.bw # timePoint = e55/e79F/e79O - only rep1 for e55
/gpfs/gibbs/pi/noonan/ak2267/NoonanBrain/mm9/ac/${timePoint}_ac_rep{1/2}.bw # timePoint = e11/e14/17F/17O
/gpfs/gibbs/pi/noonan/ak2267/NoonanBrain/mm9/me2/${timePoint}_me2_rep{1/2}.bw # timePoint = e11/e14/17F/17O
# Noonan LIMB ChIP peaks and signal:
/home/ak2267/project/EnhancerClasses/Limb/hg19/merge_${timePoint}_overlap_named.bed # timePoint = E33/E41/E44/E47
/home/ak2267/project/EnhancerClasses/Limb/rheMac2/merge_${timePoint}_overlap_named.bed # timePoint = e31/e36
/home/ak2267/project/EnhancerClasses/Limb/mm9/merge_${timePoint}_overlap_named.bed # timePoint = e10.5/e11.5/e12.5/e13.5
/gpfs/gibbs/pi/noonan/ak2267/NoonanLimb/hg19/${timePoint}_ac_rep{1/2}.bw # timePoint = E33/E41/E44/E47
/gpfs/gibbs/pi/noonan/ak2267/NoonanLimb/rheMac2/${timePoint}_ac_rep{1/2}.bw # timePoint = e31/e32/e33/e36
/gpfs/gibbs/pi/noonan/ak2267/NoonanLimb/mm9/${timePoint}_ac_rep{1/2/3}.bw # timePoint = e10.5/e11.5/e12.5/e13.5 - only rep3 for e10.5

# Use 3 prepFiles - for 1) Roller histone data, 2) Liver TF data including CTCF, 3) Noonan histone data
cd /gpfs/gibbs/pi/noonan/ak2267/speciesPairs/CGIcentric

# 1) Roller histone data
rm Roller_summaryFiles/*
python prepFiles_speciesPairs_CGIcentric_Roller.py /gpfs/gibbs/pi/noonan/ak2267/Roller/peaks/repsToUseRoller.txt repeatFilter > 230107_speciesPairs_CGIcentric_Roller_jobFile.txt
dsq --job-file 230107_speciesPairs_CGIcentric_Roller_jobFile.txt --mem-per-cpu 5G -c 1 --mail-type FAIL,END --max-jobs 30
sbatch dsq-230107_speciesPairs_CGIcentric_Roller_jobFile-2023-01-24.sh # 20240113

cd /gpfs/gibbs/pi/noonan/ak2267/speciesPairs/CGIcentric
tar -zcvf Roller_summaryFiles.gz Roller_summaryFiles/
cd /Users/acadiak/Desktop/CGI/speciesPairs/CGIcentric
scp ak2267@ruddle.hpc.yale.edu:/gpfs/gibbs/pi/noonan/ak2267/speciesPairs/CGIcentric/Roller_summaryFiles.gz .


# 2) Liver TF data including CTCF
rm LiverTF_summaryFiles/*
python prepFiles_speciesPairs_CGIcentric_LiverTF.py /gpfs/gibbs/pi/noonan/ak2267/LiverTF/peaks/repsToUseLiverTF.txt > 230109_speciesPairs_CGIcentric_LiverTF_jobFile.txt
dsq --job-file 230109_speciesPairs_CGIcentric_LiverTF_jobFile.txt --mem-per-cpu 5G -c 1 --mail-type FAIL,END --max-jobs 10
sbatch dsq-230109_speciesPairs_CGIcentric_LiverTF_jobFile-2023-01-24.sh # 20240175

cd /gpfs/gibbs/pi/noonan/ak2267/speciesPairs/CGIcentric
tar -zcvf LiverTF_summaryFiles.gz LiverTF_summaryFiles/
cd /Users/acadiak/Desktop/CGI/speciesPairs/CGIcentric
scp ak2267@ruddle.hpc.yale.edu:/gpfs/gibbs/pi/noonan/ak2267/speciesPairs/CGIcentric/LiverTF_summaryFiles.gz .


# 3) Noonan histone data
rm Noonan_summaryFiles/*
python prepFiles_speciesPairs_CGIcentric_Noonan.py repeatFilter > 230107_speciesPairs_CGIcentric_Noonan_jobFile.txt
dsq --job-file 230107_speciesPairs_CGIcentric_Noonan_jobFile.txt --mem-per-cpu 5G -c 1 --mail-type FAIL,END --max-jobs 10
sbatch dsq-230107_speciesPairs_CGIcentric_Noonan_jobFile-2023-01-24.sh # 20240194

cd /gpfs/gibbs/pi/noonan/ak2267/speciesPairs/CGIcentric
tar -zcvf Noonan_summaryFiles.gz Noonan_summaryFiles/
cd /Users/acadiak/Desktop/CGI/speciesPairs/CGIcentric
scp ak2267@ruddle.hpc.yale.edu:/gpfs/gibbs/pi/noonan/ak2267/speciesPairs/CGIcentric/Noonan_summaryFiles.gz .


######## FOR FIGURE 2
# Run all species pairs with H3K27ac, but remove step that filters oCGIs overlapping peaks associated with promoters and other features
# Keep filters based on length and repeat content
# The resulting tables will be used for Figure 2 to assess oCGI numbers across species (the bar plots in Fig 2B)
# And to assess differences like length and CpG number

cd /gpfs/gibbs/pi/noonan/ak2267/speciesPairs/CGIcentric

mkdir SummaryFiles_noPromFilter # write output to here
rm SummaryFiles_noPromFilter/*

python prepFiles_speciesPairs_CGIcentric_Roller_noPromFilter.py /gpfs/gibbs/pi/noonan/ak2267/Roller/peaks/repsToUseRoller.txt repeatFilter > 230110_speciesPairs_CGIcentric_Roller_noPromFilter_jobFile.txt
dsq --job-file 230110_speciesPairs_CGIcentric_Roller_noPromFilter_jobFile.txt --mem-per-cpu 5G -c 1 --mail-type FAIL,END --max-jobs 10
sbatch dsq-230110_speciesPairs_CGIcentric_Roller_noPromFilter_jobFile-2023-01-24.sh # 20240243

python prepFiles_speciesPairs_CGIcentric_Noonan_noPromFilter.py repeatFilter > 230110_speciesPairs_CGIcentric_Noonan_noPromFilter_jobFile.txt
dsq --job-file 230110_speciesPairs_CGIcentric_Noonan_noPromFilter_jobFile.txt --mem-per-cpu 5G -c 1 --mail-type FAIL,END
sbatch dsq-230110_speciesPairs_CGIcentric_Noonan_noPromFilter_jobFile-2023-01-24.sh # 20240259

# download for use in R (Fig2_rhesusMouse.R)
cd /gpfs/gibbs/pi/noonan/ak2267/speciesPairs/CGIcentric
tar -zcvf SummaryFiles_noPromFilter.gz SummaryFiles_noPromFilter/
cd /Users/acadiak/Desktop/CGI/speciesPairs/CGIcentric
scp ak2267@ruddle.hpc.yale.edu:/gpfs/gibbs/pi/noonan/ak2267/speciesPairs/CGIcentric/SummaryFiles_noPromFilter.gz .


##### put reconciledCGIs and reconciledPeaks on lab server for viewing on the browser
# 8/12/22

ssh ak2267@10.5.37.220
cd /home/ak2267/akocher_www/CGI
mkdir hg19
mkdir rheMac2
mkdir mm9

cd /gpfs/gibbs/pi/noonan/ak2267/speciesPairs/CGIcentric

# Roller data first - turn reconciledCGI and reconciledPeak files for each species pair into bigBeds, and transfer to appropriate folder on server
# will need to be incorporated in script that makes trackDb.txt files ()

mkdir bigBed
for speciesA in rheMac10 calJac4 mm39 rn7 susScr11 canFam6 felCat9 equCab3; do for speciesB in calJac4 mm39 rn7 susScr11 canFam6 felCat9 equCab3; do for tissue in brain liver muscle testis; do for mark in H3K4me3 H3K27ac H3K4me1; do dir=Roller_${speciesA}_${speciesB}_${tissue}_${mark}; if [ -d "$dir" ]; then echo $dir ; sort -k1,1 -k2,2n ${dir}/speciesA_reconciledCGI_noFeaturePeaks.bed > bigBed/temp_sort.bed ; bedToBigBed bigBed/temp_sort.bed /home/ak2267/genomes/chrom.sizes/${speciesA}.chrom.sizes bigBed/reconciledCGIs_${speciesA}_with${speciesB^}_${tissue}_${mark}.bb ; sort -k1,1 -k2,2n ${dir}/speciesB_reconciledCGI_noFeaturePeaks.bed > bigBed/temp_sort.bed ; bedToBigBed bigBed/temp_sort.bed /home/ak2267/genomes/chrom.sizes/${speciesB}.chrom.sizes bigBed/reconciledCGIs_${speciesB}_with${speciesA^}_${tissue}_${mark}.bb ; scp bigBed/reconciledCGIs_${speciesA}_with${speciesB^}_${tissue}_${mark}.bb ak2267@10.5.37.220:/home/ak2267/akocher_www/CGI/${speciesA} ; scp bigBed/reconciledCGIs_${speciesB}_with${speciesA^}_${tissue}_${mark}.bb ak2267@10.5.37.220:/home/ak2267/akocher_www/CGI/${speciesB} ; sort -k1,1 -k2,2n ${dir}/speciesA_reconciledPeaks.bed > bigBed/temp_sort.bed; bedToBigBed bigBed/temp_sort.bed /home/ak2267/genomes/chrom.sizes/${speciesA}.chrom.sizes bigBed/reconciledPeaks_${speciesA}_with${speciesB^}_${tissue}_${mark}.bb ; sort -k1,1 -k2,2n ${dir}/speciesB_reconciledPeaks.bed > bigBed/temp_sort.bed ; bedToBigBed bigBed/temp_sort.bed /home/ak2267/genomes/chrom.sizes/${speciesB}.chrom.sizes bigBed/reconciledPeaks_${speciesB}_with${speciesA^}_${tissue}_${mark}.bb ; scp bigBed/reconciledPeaks_${speciesA}_with${speciesB^}_${tissue}_${mark}.bb ak2267@10.5.37.220:/home/ak2267/akocher_www/CGI/${speciesA} ; scp bigBed/reconciledPeaks_${speciesB}_with${speciesA^}_${tissue}_${mark}.bb ak2267@10.5.37.220:/home/ak2267/akocher_www/CGI/${speciesB}; sleep 20s; fi;done;done;done;done

# modify makeTrackDb.py to include these new files for Roller data - default them to be hidden
cd /gpfs/gibbs/pi/noonan/ak2267/trackDb
python makeTrackDb.py
for species in rheMac10 calJac4 mm39 rn7 susScr11 canFam6 felCat9 equCab3
do
    scp trackDb_${species}.txt ak2267@10.5.37.220:/home/ak2267/akocher_www/CGI/${species}
done

# add hg19, rheMac2, and mm9 to genomes.txt
cd /home/ak2267/project/Roller/ChIP/batch
scp genomes.txt ak2267@10.5.37.220:/home/ak2267/akocher_www/CGI

# move Noonan reconciledCGIs and reconciledPeaks to the browser
cd /gpfs/gibbs/pi/noonan/ak2267/speciesPairs/CGIcentric
for speciesA in hg19 rheMac2; do for speciesB in rheMac2 mm9; do for tissue in brain limb; do for mark in ac me2; do for timePoint in 0 1 2 3; do dir=Noonan_${speciesA}_${speciesB}_${tissue}_${mark}_${timePoint}; if [ -d "$dir" ]; then echo $dir ; sort -k1,1 -k2,2n ${dir}/speciesA_reconciledCGI_noFeaturePeaks.bed > bigBed/temp_sort.bed ; bedToBigBed bigBed/temp_sort.bed /home/ak2267/genomes/chrom.sizes/${speciesA}.chrom.sizes bigBed/reconciledCGIs_${speciesA}_with${speciesB^}_${tissue}_${mark}_${timePoint}.bb ; sort -k1,1 -k2,2n ${dir}/speciesB_reconciledCGI_noFeaturePeaks.bed > bigBed/temp_sort.bed ; bedToBigBed bigBed/temp_sort.bed /home/ak2267/genomes/chrom.sizes/${speciesB}.chrom.sizes bigBed/reconciledCGIs_${speciesB}_with${speciesA^}_${tissue}_${mark}_${timePoint}.bb ; scp bigBed/reconciledCGIs_${speciesA}_with${speciesB^}_${tissue}_${mark}_${timePoint}.bb ak2267@10.5.37.220:/home/ak2267/akocher_www/CGI/${speciesA} ; scp bigBed/reconciledCGIs_${speciesB}_with${speciesA^}_${tissue}_${mark}_${timePoint}.bb ak2267@10.5.37.220:/home/ak2267/akocher_www/CGI/${speciesB};sort -k1,1 -k2,2n ${dir}/speciesA_reconciledPeaks.bed > bigBed/temp_sort.bed; bedToBigBed bigBed/temp_sort.bed /home/ak2267/genomes/chrom.sizes/${speciesA}.chrom.sizes bigBed/reconciledPeaks_${speciesA}_with${speciesB^}_${tissue}_${mark}_${timePoint}.bb ; sort -k1,1 -k2,2n ${dir}/speciesB_reconciledPeaks.bed > bigBed/temp_sort.bed ; bedToBigBed bigBed/temp_sort.bed /home/ak2267/genomes/chrom.sizes/${speciesB}.chrom.sizes bigBed/reconciledPeaks_${speciesB}_with${speciesA^}_${tissue}_${mark}_${timePoint}.bb ; scp bigBed/reconciledPeaks_${speciesA}_with${speciesB^}_${tissue}_${mark}_${timePoint}.bb ak2267@10.5.37.220:/home/ak2267/akocher_www/CGI/${speciesA} ; scp bigBed/reconciledPeaks_${speciesB}_with${speciesA^}_${tissue}_${mark}_${timePoint}.bb ak2267@10.5.37.220:/home/ak2267/akocher_www/CGI/${speciesB};sleep 20s;fi;done;done;done;done;done

# move Noonan UCSC AL CGIs to the server
cd /home/ak2267/genomes/CGI/UCSC_AL
for species in hg19 rheMac2 mm9
do
    cut -f 1,2,3 ${species}_CGIsAL.bed > ${species}_3col.bed ; bedToBigBed ${species}_3col.bed /home/ak2267/genomes/chrom.sizes/${species}.chrom.sizes ${species}_CGIsAL.bb
    scp ${species}_CGIsAL.bb ak2267@10.5.37.220:/home/ak2267/akocher_www/CGI/${species}
done

# move human phastCons elements bed file to the server
cd /home/ak2267/genomes/phastCons
sort -k1,1 -k2,2n phastConsElements100way_hg19.bed > phastConsElements100way_hg19_sorted.bed
bedToBigBed phastConsElements100way_hg19_sorted.bed /home/ak2267/genomes/chrom.sizes/hg19.chrom.sizes phastConsElements100way_hg19.bb
scp phastConsElements100way_hg19.bb ak2267@10.5.37.220:/home/ak2267/akocher_www/CGI/hg19

# run python script to make trackDb.txt files and transfer to the server
cd /gpfs/gibbs/pi/noonan/ak2267/trackDb
python makeTrackDb_Noonan.py
for species in hg19 rheMac2 mm9
do
    scp trackDb_${species}.txt ak2267@10.5.37.220:/home/ak2267/akocher_www/CGI/${species}
done

