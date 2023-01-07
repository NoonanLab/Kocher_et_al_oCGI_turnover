Pipeline for analysis of oCGIs and histone modification peaks between species pairs

#### oCGI-CENTRIC ANALYSIS ####

220726_speciesPairs.sh
# Generate a table for each species pair x tissue x mark x time point for analysis downstream in R
# Uses several python scripts to make slurm job files:
prepFiles_speciesPairs_CGIcentric_Roller.py
prepFiles_speciesPairs_CGIcentric_LiverTF.py # ****need to find where I made repsToUseLiverTF.txt****
prepFiles_speciesPairs_CGIcentric_Noonan.py

221017_runPermutationOnHPC.sh
# Run permutation test to assess enrichment in each box in the grid using R on the cluster

#### PEAK-CENTRIC ANALYSIS ####

221130_speciesPairs_peakCentric.sh
# Generate a table for each species pair x tissue x mark x time point for analysis downstream in R
# But this time the analysis is PEAK-CENTRIC instead of oCGI-centric
# Uses several python scripts to make slurm job files:
prepFiles_speciesPairs_peakCentric_Roller.py
prepFiles_speciesPairs_peakCentric_Noonan.py

230109_runPermutationOnHPC_peakCentric.sh
# Run permutation test to assess enrichment in each box in the grid using R on the cluster
# For the PEAK-CENTRIC pipeline
