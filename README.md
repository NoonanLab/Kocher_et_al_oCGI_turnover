# Analysis pipeline for Kocher et al (2023)
CpG island turnover events predict evolutionary changes in enhancer activity  
bioRxiv https://doi.org/10.1101/2023.05.09.540063  

This repository contains a separate directory for each step of the analysis  
&emsp;&emsp;Step 1 - identify oCGIs in nine mammalian genomes  
&emsp;&emsp;Step 2 - process published ChIP-seq data  
&emsp;&emsp;Step 3 - overlap of oCGIs and ChIP-seq peaks in single species  
&emsp;&emsp;Step 4 - identify orthologous oCGIs within species pairs  
&emsp;&emsp;Step 5 - analyze oCGI turnover and peak turnover across species pairs  
&emsp;&emsp;Step 6 - enrichment of oCGI species patterns in Human Gain Enhancers (HGEs)  
&emsp;&emsp;Step 7 - process ChIP-seq and RNA-seq from the hs754 humanized mouse  
&emsp;&emsp;Step 8 - gene expression changes associated with oCGI and peak turnover  
&emsp;&emsp;Step 9 - GC-biased gene conversion  
&emsp;&emsp;Step 10 - analysis of VISTA element LacZ assay results  
&emsp;&emsp;Step 11 - generate figures and tables in R  

Each directory contains:  
&emsp;&emsp;README.md file explaining its contents  
&emsp;&emsp;Workflow files (ending .sh) for manipulating functional genomic files and running jobs from the command line  
&emsp;&emsp;Python scripts (ending .py) and R scripts (ending .R) required for steps of the analysis, called within .sh files  
&emsp;&emsp;R scripts that generate figures and tables (see Step11)
