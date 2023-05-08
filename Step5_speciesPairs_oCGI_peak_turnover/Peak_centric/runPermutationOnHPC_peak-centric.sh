# 1/9/23
# Purpose: run permutation test (20,000 permutations) on Ruddle
# with a separate job for each combination (speciesPair x tissue x mark x timePoint)
# this script is for PEAK-CENTRIC analysis (not CGI-centric)

# R code: speciesPairs_permutation_HPC_peakCentric.R
cd /gpfs/gibbs/pi/noonan/ak2267/speciesPairs/peakCentric/permutation

# need to input file path to summaryFile for each combination
# and label from which to take 

# make job file for Roller and Noonan data
for summaryFile in $(ls /gpfs/gibbs/pi/noonan/ak2267/speciesPairs/peakCentric/*_summaryFiles/*.txt)
do
    echo 'cd /gpfs/gibbs/pi/noonan/ak2267/speciesPairs/peakCentric/permutation ; source ~/.bashrc ; source ~/.bash_profile ; module load R ; Rscript speciesPairs_permutation_HPC_peak-centric.R '${summaryFile}
done >> 230109_permutation_peakCentric_jobFile.txt

# run jobs
dsq --job-file 230109_permutation_peakCentric_jobFile.txt --mem-per-cpu 5G -c 1 --mail-type FAIL,END
sbatch dsq-230109_permutation_peakCentric_jobFile-2023-01-09.sh # 19832421
sbatch dsq-230109_permutation_peakCentric_jobFile-2023-03-08.sh # 21674933

# cat results into single file for further use in R
cat *_3x3.txt > permutationResults_3x3.txt
cat *_4x3.txt > permutationResults_4x3.txt

# download to laptop for further use in R Studio
cd /Users/acadiak/Desktop/CGI/speciesPairs/peakCentric/permutation
scp ak2267@ruddle.hpc.yale.edu:/gpfs/gibbs/pi/noonan/ak2267/speciesPairs/peakCentric/permutation/permutationResults_3x3.txt .
scp ak2267@ruddle.hpc.yale.edu:/gpfs/gibbs/pi/noonan/ak2267/speciesPairs/peakCentric/permutation/permutationResults_4x3.txt .

# then see R code in qValues_peakCentric.R to generate q values (separately for 3x3 and 4x3)
# final results are in /Users/acadiak/Desktop/CGI/speciesPairs/CGIcentric/permutation
# permutation_peakCentric_qAndFD_3x3.txt
# permutation_peakCentric_qAndFD_4x3.txt