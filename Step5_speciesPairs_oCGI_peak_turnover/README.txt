Pipeline for analysis of oCGIs and histone modification peaks between species pairs

#### oCGI-CENTRIC ANALYSIS ####

220726_speciesPairs.sh
# Generate a table for each species pair x tissue x mark x time point for analysis downstream in R
# Uses several python scripts to make slurm job files:
prepFiles_speciesPairs_CGIcentric_Roller.py
prepFiles_speciesPairs_CGIcentric_LiverTF.py
prepFiles_speciesPairs_CGIcentric_Noonan.py
# All of these job files create intermediate files, then call a python script that integrates them into a final table:
countFinalOverlaps_speciesPairs_CGIcentric.py

# Make files for trackHub with all ChIP-seq data, plus all CGIs and orthologous oCGIs
makeTrackDb.py
makeTrackDb_Noonan.py

221017_runPermutationOnHPC.sh
# Run permutation test to assess enrichment in each box in the grid using R on the cluster

qValues_CGIcentric.R
# Analyze permutation results in R - convert p values to q values with correction for multiple testing, and report Fold Differences (FD)

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
