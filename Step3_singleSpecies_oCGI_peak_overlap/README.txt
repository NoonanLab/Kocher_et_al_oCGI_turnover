
# Analysis in single species for Fig 1
220912_CGI_singleSpecies_activity.sh

######### Make summary tables for CGI-centric analysis
# Commands are written out within 220912_CGI_singleSpecies_activity.sh

# But pipeline calls this python script to make commands for featureCounts
makeCommands_singleSpeciesFeatureCounts.py

# Also calls this python script to summarize all intermediate files at end of pipeline
makeSingleSpeciesSummaryTables_CGIcentric.py

######### Make summary tables for Peak-centric analysis
# Pipeline has been moved from .sh file to job files made by this python script
prepFiles_singleSpecies_peakCentric.py

# Which calls this python script to summarize all intermediate files at end of pipeline
makeSummary_singleSpecies_peakCentric.py

# Also calls the following python scripts that are in general python_scripts directory
makeGTF_nameFromCol4.py
countRPM_general.py
countBigWigSignal_general.py

######### Determine enrichment by shuffling sites on the genome background
# Use this prepFiles python script to make job files
prepFiles_singleSpecies_shuffling.py

# Use these python scripts to summarize the results
summarizeReshuffling_CGI-centric.py
summarizeReshuffling_peak-centric.py
