# Step 8: RNA changes associated with oCGI turnover

### Analyze RNA-seq data from Roller, Stamper et al 2021 and perform test of association with oCGI turnover
RollerRNA.sh

### Calls the following:  
convertChrNamesToUCSCinGTF.py  
&emsp;&emsp;convert Ensembl chromosome names to be compatible with UCSC (add "chr" before numbers)  
prepFiles_RollerRNA_Nov22.py  
&emsp;&emsp;process RNA-seq data: download, align, count reads in genes  
outputCGIsAsBed.py  
&emsp;&emsp;within each species pair, tissue, and mark:  
&emsp;&emsp;make bed files with (1) A-only oCGIs in A-only peaks, (2) with B-only oCGIs in B-only peaks  
&emsp;&emsp;where names in column 4 tell their species-specificity  
makeCategoryBedFiles_ENSEMBL.py  
&emsp;&emsp;input GTF file, output bed files for all feature types; want to use protein-coding promoters  
extendRegulatoryDomains.py  
&emsp;&emsp;turn promoter bed file into bed file with regulatory domains, as defined by GREAT rules  
makeTwoWayOrthologTable.py  
&emsp;&emsp;take all info and summarize into a table for each species pair, tissue, and mark  

## Perform resampling test using the tables made by makeTwoWayOrthologTable.py
RNA_resamplingTest.R  
&emsp;&emsp;perform resampling test  
summarizeRNAresampling.py  
&emsp;&emsp;summarize results of test for further use in R
