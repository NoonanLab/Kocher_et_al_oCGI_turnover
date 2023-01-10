# 11/30/22
# Purpose: model this pipeline after 220726_speciesPairs.sh
# but invert so that analysis is peak-centric instead of CGI-centric
# Roller species: H3K4me3, H3K27ac, H3K4me1, CTCF, 4 Liver TFs
# Noonan species: Reilly H3K27ac & H3K4me2, Cotney H3K27ac

# Reconciled CGIs are here:
/gpfs/gibbs/pi/noonan/ak2267/speciesPairs/consensusCGIs/{speciesA}_{speciesB}/{speciesA}_{speciesB}_speciesA_reconciledCGI.bed
/gpfs/gibbs/pi/noonan/ak2267/speciesPairs/consensusCGIs/{speciesA}_{speciesB}/{speciesA}_{speciesB}_speciesB_reconciledCGI.bed

################ Peak-centric
# Work here:
mkdir /gpfs/gibbs/pi/noonan/ak2267/speciesPairs/peakCentric
cd /gpfs/gibbs/pi/noonan/ak2267/speciesPairs/peakCentric

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

# Use 2 prepFiles - for 1) Roller histone data, 2) Noonan histone data
cd /gpfs/gibbs/pi/noonan/ak2267/speciesPairs/peakCentric

# 1) Roller histone data
rm Roller_summaryFiles/*
python prepFiles_speciesPairs_peakCentric_Roller.py /gpfs/gibbs/pi/noonan/ak2267/Roller/peaks/repsToUseRoller.txt noRepeatFilter > 221201_speciesPairs_peakCentric_Roller_jobFile.txt
dsq --job-file 221201_speciesPairs_peakCentric_Roller_jobFile.txt --mem-per-cpu 5G -c 1 --mail-type FAIL,END
sbatch dsq-221201_speciesPairs_peakCentric_Roller_jobFile-2022-12-02.sh # 18878176

# 2) Noonan histone data
python prepFiles_speciesPairs_peakCentric_Noonan.py noRepeatFilter > 221202_speciesPairs_peakCentric_Noonan_jobFile.txt
dsq --job-file 221202_speciesPairs_peakCentric_Noonan_jobFile.txt --mem-per-cpu 5G -c 1 --mail-type FAIL,END
sbatch dsq-221202_speciesPairs_peakCentric_Noonan_jobFile-2022-12-02.sh # 18880962

# results are here
/gpfs/gibbs/pi/noonan/ak2267/speciesPairs/peakCentric/Roller_summaryFiles
/gpfs/gibbs/pi/noonan/ak2267/speciesPairs/peakCentric/Noonan_summaryFiles

# RUN PERMUTATION TEST WITH R ON RUDDLE
# see 221017_runPermutationOnHPC_peakCentric.sh
# which calls speciesPairs_permutation_HPC_peakCentric.R
# then use qValues_peakCentric.R to make permutation_peakCentric_qAndFD_3x3.txt and permutation_peakCentric_qAndFD_4x3.txt for use in R for plots
