# Step 3: overlap of oCGIs and peaks within single species

## Make summary tables for CGI-centric analysis (for use in R)
singleSpecies_oCGIs_peaks.sh  
&emsp;&emsp;Calls makeCommands_singleSpeciesFeatureCounts.py  
&emsp;&emsp;Calls makeSingleSpeciesSummaryTables_CGIcentric.py

## Determine level of oCGI enrichment for peaks by shuffling sites on the genome background
singleSpecies_oCGIs_peaks.sh  
&emsp;&emsp;Calls prepFiles_singleSpecies_shuffling.py  
&emsp;&emsp;Calls summarizeReshuffling_CGI-centric.py

## Make summary tables for Peak-centric analysis (for use in R)
singleSpecies_oCGIs_peaks.sh  
&emsp;&emsp;Calls prepFiles_singleSpecies_peakCentric.py  
&emsp;&emsp;Calls makeSummary_singleSpecies_peakCentric.py

## These scripts also reference the following that are in the python_scripts directory in the main repository:
makeGTF_nameFromCol4.py  
countRPM_general.py  
countBigWigSignal_general.py  
restrictToSitesThatMapBack.py  
restrictToLO.py  
