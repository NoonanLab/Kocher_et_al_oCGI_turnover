# 3/7/2024
# Purpose: re-run pipeline using a difference threshold for species-specific CGIs

# I have to re-run my pipeline, adding a step where I capture G and C counts in each reconciled CGI
# The existing pipeline reports GC percent, but not individual counts of C and G
# which I need to calculate Obs/Exp ratios for both species

# I will need a new version of countFinalOverlaps_speciesPairs_CGIcentric.py
# and to run it using commands made by prepFiles_speciesPairs_CGIcentric_Roller_noPromFilter.py

# new versions of these files:
# ak20240307_countFinalOverlaps_speciesPairs_CGIcentric_addGCcounts.py
# ak20240307_prepFiles_speciesPairs_CGIcentric_Roller_noPromFilter_addGCcounts.py

# commands are modified from speciesPairs_oCGI-centric.sh


######## USE SAME PIPELINE AS FOR FIGURE 2
# Run all species pairs with H3K27ac, but remove step that filters oCGIs overlapping peaks associated with promoters and other features
# Keep filters based on length and repeat content
# The resulting tables will be used for Figure 2 to assess oCGI numbers across species (the bar plots in Fig 2B)
# And to assess differences like length and CpG number

cd /gpfs/gibbs/pi/noonan/ak2267/speciesPairs/CGIcentric

mkdir SummaryFiles_noPromFilter_addGCcounts # write output to here
rm SummaryFiles_noPromFilter_addGCcounts/*

python ak20240307_prepFiles_speciesPairs_CGIcentric_Roller_noPromFilter_addGCcounts.py /gpfs/gibbs/pi/noonan/ak2267/Roller/peaks/repsToUseRoller.txt repeatFilter > 240307_speciesPairs_CGIcentric_Roller_noPromFilter_jobFile.txt
dsq --job-file 240307_speciesPairs_CGIcentric_Roller_noPromFilter_jobFile.txt --mem-per-cpu 5G -c 1 --mail-type FAIL,END --max-jobs 10
sbatch dsq-240307_speciesPairs_CGIcentric_Roller_noPromFilter_jobFile-2024-03-07.sh # 23031875
# ^^ this made files with messed up headers, need to re-run the final python step with an edited version of the countFinalOverlaps script

rm SummaryFiles_noPromFilter_addGCcounts/*

python ak20240307_pythonQuickFix.py /gpfs/gibbs/pi/noonan/ak2267/Roller/peaks/repsToUseRoller.txt repeatFilter > 240307_pythonStepOnly_jobFile.txt
dsq --job-file 240307_pythonStepOnly_jobFile.txt --mem-per-cpu 5G -c 1 --mail-type FAIL,END --max-jobs 10
sbatch dsq-240307_pythonStepOnly_jobFile-2024-03-07.sh # 23040007

# Don't do Noonan data, only Roller

# download for use in R (ak20240307_differenceRequirement.R)
cd /gpfs/gibbs/pi/noonan/ak2267/speciesPairs/CGIcentric
tar -zcvf SummaryFiles_noPromFilter_addGCcounts.gz SummaryFiles_noPromFilter_addGCcounts/
cd /Users/acadiak/Desktop/Yale/\!Writing/Kocher_CpG/GB_Revisions/Revised_code/Step12_revisions/Difference_requirement
scp ak2267@mccleary.ycrc.yale.edu:/gpfs/gibbs/pi/noonan/ak2267/speciesPairs/CGIcentric/SummaryFiles_noPromFilter_addGCcounts.gz .




######## RUN AGAIN BUT THIS TIME USE FACOUNT TO GET COUNTS IN ORIGINAL CGIs, NOT RECONCILED
# for use in addressing reviewer question about Fig 2C and S12

# make new prepFiles: ak20240307_prepFiles_speciesPairs_CGIcentric_Roller_noPromFilter_OriginalGCcounts.py

cd /gpfs/gibbs/pi/noonan/ak2267/speciesPairs/CGIcentric

mkdir SummaryFiles_noPromFilter_OriginalGCcounts # write output to here
rm SummaryFiles_noPromFilter_OriginalGCcounts/*

python ak20240307_prepFiles_speciesPairs_CGIcentric_Roller_noPromFilter_OriginalGCcounts.py /gpfs/gibbs/pi/noonan/ak2267/Roller/peaks/repsToUseRoller.txt repeatFilter > 240307_speciesPairs_CGIcentric_Roller_noPromFilter_OriginalGCcounts_jobFile.txt
dsq --job-file 240307_speciesPairs_CGIcentric_Roller_noPromFilter_OriginalGCcounts_jobFile.txt --mem-per-cpu 5G -c 1 --mail-type FAIL,END --max-jobs 10
sbatch dsq-240307_speciesPairs_CGIcentric_Roller_noPromFilter_OriginalGCcounts_jobFile-2024-03-07.sh # 23094309
# run again 20240321 but with fasta from MASKED genomes
rm SummaryFiles_noPromFilter_OriginalGCcounts/*
python ak20240307_prepFiles_speciesPairs_CGIcentric_Roller_noPromFilter_OriginalGCcounts.py /gpfs/gibbs/pi/noonan/ak2267/Roller/peaks/repsToUseRoller.txt repeatFilter > 240307_speciesPairs_CGIcentric_Roller_noPromFilter_OriginalGCcounts_jobFile.txt
dsq --job-file 240307_speciesPairs_CGIcentric_Roller_noPromFilter_OriginalGCcounts_jobFile.txt --mem-per-cpu 5G -c 1 --mail-type FAIL,END
sbatch dsq-240307_speciesPairs_CGIcentric_Roller_noPromFilter_OriginalGCcounts_jobFile-2024-03-21.sh # 24827652


# download and plot in ak20240307_Fig2_Obs_Exp.R
cd /gpfs/gibbs/pi/noonan/ak2267/speciesPairs/CGIcentric
tar -zcvf SummaryFiles_noPromFilter_OriginalGCcounts.gz SummaryFiles_noPromFilter_OriginalGCcounts/
cd /Users/acadiak/Desktop/Yale/\!Writing/Kocher_CpG/GB_Revisions/Revised_code/Step12_revisions/Difference_requirement
scp ak2267@mccleary.ycrc.yale.edu:/gpfs/gibbs/pi/noonan/ak2267/speciesPairs/CGIcentric/SummaryFiles_noPromFilter_OriginalGCcounts.gz .




### RUN PIPELINE FOR FIGURES 3 AND 4 GRIDS USING MODIFIED SCRIPTS TO RESTRICT BASED ON DIFFERENCE REQUIREMENT
# use 3 different cutoffs: 0.1, 0.3, and 0.5
# generate for all species although I only intend to make plots for the subset that I show above
# the cutoffs will be applied in R in the permutation tests I run on the cluster
# what I need to re-run here first is generating the tables for R using the modifications I made to countFinalOverlap script above
# so that it captures G and C counts, so that in R I can calculate Obs/Exp ratios

cd /gpfs/gibbs/pi/noonan/ak2267/speciesPairs/CGIcentric

mkdir Roller_summaryFiles_addGCcounts/
rm Roller_summaryFiles_addGCcounts/*

# need a new prepFiles: ak20240307_prepFiles_speciesPairs_CGIcentric_Roller_addGCcounts.py
# modified to add back in filter based on overlapping a promoter peak - reference original version, prepFiles_speciesPairs_CGIcentric_Roller.py
# that calls the countFinalOverlaps I made above: ak20240307_countFinalOverlaps_speciesPairs_CGIcentric_addGCcounts.py

python ak20240307_prepFiles_speciesPairs_CGIcentric_Roller_addGCcounts.py /gpfs/gibbs/pi/noonan/ak2267/Roller/peaks/repsToUseRoller.txt repeatFilter > 240307_speciesPairs_CGIcentric_Roller_jobFile.txt
dsq --job-file 240307_speciesPairs_CGIcentric_Roller_jobFile.txt --mem-per-cpu 5G -c 1 --mail-type FAIL,END --max-jobs 100 # 
sbatch dsq-240307_speciesPairs_CGIcentric_Roller_jobFile-2024-03-07.sh # 23072995

# if you wanted to download these:
#cd /gpfs/gibbs/pi/noonan/ak2267/speciesPairs/CGIcentric
#tar -zcvf Roller_summaryFiles_addGCcounts.gz Roller_summaryFiles_addGCcounts/
#cd /Users/acadiak/Desktop/CGI/speciesPairs/CGIcentric
#scp ak2267@mccleary.ycrc.hpc.yale.edu:/gpfs/gibbs/pi/noonan/ak2267/speciesPairs/CGIcentric/Roller_summaryFiles_addGCcounts.gz .


# run permutation on HPC using each difference threshold - model after runPermutation on HPC_oCGI-centric.sh
# which calls speciesPairs_permutation_HPC_oCGI-centric.R
# make new version: ak20240307_speciesPairs_permutation_HPC_oCGI-centric_differenceRequirement.R

mkdir permutation_diffReq

for summaryFile in $(ls /gpfs/gibbs/pi/noonan/ak2267/speciesPairs/CGIcentric/Roller_summaryFiles_addGCcounts/*.txt)
do
    for diffReq in '0' '0.1' '0.3' '0.5'
    do
        echo 'cd /gpfs/gibbs/pi/noonan/ak2267/speciesPairs/CGIcentric/permutation_diffReq ; source ~/.bashrc ; source ~/.bash_profile ; module load R ; Rscript ak20240307_speciesPairs_permutation_HPC_oCGI-centric_differenceRequirement.R '${summaryFile} ${diffReq}
    done
done >> 240307_permutation_CGIcentric_diffReq_jobFile.txt

dsq --job-file 240307_permutation_CGIcentric_diffReq_jobFile.txt --mem-per-cpu 5G -c 1 --mail-type FAIL,END --max-jobs 100
sbatch dsq-240307_permutation_CGIcentric_diffReq_jobFile-2024-03-07.sh # 23074657

# note that rheMac10_calJac4 and rheMac10_equCab3 failed for H3K4me3 in the liver - possibly due to not enough sites? did not investigate

# cat results into single file for further use in R
cd permutation_diffReq
cat *_3x3*.txt > permutationResults_diffReq_3x3.txt

# then download
# see R code in ak20240307_qValues_CGIcentric.R to generate q values (correct across all tests done for these grids today)
# final results are in /Users/acadiak/Desktop/Yale/!Writing/Kocher_CpG/GB_Revisions/Revised_code/Step12_revisions/Difference_requirement
# permutation_qAndFD_3x3.txt

cd /Users/acadiak/Desktop/Yale/\!Writing/Kocher_CpG/GB_Revisions/Revised_code/Step12_revisions/Difference_requirement
scp ak2267@mccleary.ycrc.yale.edu:/gpfs/gibbs/pi/noonan/ak2267/speciesPairs/CGIcentric/permutation_diffReq/permutationResults_diffReq_3x3.txt .




