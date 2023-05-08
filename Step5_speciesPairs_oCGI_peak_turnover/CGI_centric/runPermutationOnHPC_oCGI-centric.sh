# 10/17/22
# Purpose: run permutation test (20,000 permutations) on Ruddle
# with a separate job for each combination (speciesPair x tissue x mark x timePoint)

# R code: speciesPairs_permutation_HPC_CGIcentric.R
cd /gpfs/gibbs/pi/noonan/ak2267/speciesPairs/CGIcentric/permutation

# need to input file path to summaryFile for each combination
# and label from which to take 

# make job file for Roller, Noonan, and LiverTF data (analyze all at once, and correct p-values in R all at once - see below)
for summaryFile in $(ls /gpfs/gibbs/pi/noonan/ak2267/speciesPairs/CGIcentric/*_summaryFiles/*.txt)
do
    echo 'cd /gpfs/gibbs/pi/noonan/ak2267/speciesPairs/CGIcentric/permutation ; source ~/.bashrc ; source ~/.bash_profile ; module load R ; Rscript speciesPairs_permutation_HPC_oCGI-centric.R '${summaryFile}
done >> 230109_permutation_CGIcentric_jobFile.txt

# run jobs
dsq --job-file 230109_permutation_CGIcentric_jobFile.txt --mem-per-cpu 5G -c 1 --mail-type FAIL,END --max-jobs 100
sbatch dsq-230109_permutation_CGIcentric_jobFile-2023-01-24.sh # 20241555
sbatch dsq-230109_permutation_CGIcentric_jobFile-2023-03-05.sh # 21637978

# cat results into single file for further use in R
cat *_3x3.txt > permutationResults_3x3.txt
cat *_4x3.txt > permutationResults_4x3.txt

# download to laptop for further use in R Studio
cd /Users/acadiak/Desktop/CGI/speciesPairs/CGIcentric/permutation
scp ak2267@ruddle.hpc.yale.edu:/gpfs/gibbs/pi/noonan/ak2267/speciesPairs/CGIcentric/permutation/permutationResults_3x3.txt .
scp ak2267@ruddle.hpc.yale.edu:/gpfs/gibbs/pi/noonan/ak2267/speciesPairs/CGIcentric/permutation/permutationResults_4x3.txt .

# then see R code in qValues_CGIcentric.R to generate q values (separately for 3x3 and 4x3)
# final results are in /Users/acadiak/Desktop/CGI/speciesPairs/CGIcentric/permutation
# permutation_qAndFD_3x3.txt
# permutation_qAndFD_4x3.txt
