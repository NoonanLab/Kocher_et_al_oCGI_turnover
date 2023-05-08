# Step 8: RNA changes associated with oCGI turnover

### Analyze RNA-seq data from Roller, Stamper et al 2021 and perform test of association with oCGI turnover
RollerRNA.sh

### Calls the following:  
convertChrNamesToUCSCinGTF.py - convert Ensembl chromosome names to be compatible with UCSC (add "chr" before numbers)  
assessOrthologDistribution.py - count 1:1 orthologs that fall on unassembled contigs - very few, so I exclude them  
assessNumberOfOrthologs.py - count 1:1 orthologs that are 1:1 instead of 1:many or many:1 or many:many - minority, so I exclude them  
prepFiles_RollerRNA_Nov22.py - process RNA-seq data: download, align, count reads in genes  
outputCGIsAsBed.py - within each species pair, tissue, and mark:  
&emsp;&emsp;make bed file with A-only oCGIs in A-only peaks, and a bed file with B-only oCGIs in B-only peaks, with names in column 4 telling their species-specificity  
outputUniqueBiotype.py - output a list of all unique biotypes in Ensembl GTF files, for identifying protein-coding genes in makeCategoryBedFiles_ENSEMBL.py  
makeCategoryBedFiles_ENSEMBL.py - input GTF file, output bed files for all feature types; want to use protein-coding promoters  
extendRegulatoryDomains.py - turn promoter bed file into bed file with regulatory domains, as defined by GREAT rules  
makeTwoWayOrthologTable.py - take all info and summarize into a table for each species pair, tissue, and mark  

## Perform resampling test using the tables made by makeTwoWayOrthologTable.py
RNA_resamplingTest.R - perform resampling test  
summarizeRNAresampling.py - summarize results of test for further use in R
