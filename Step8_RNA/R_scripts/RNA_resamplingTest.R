# 2/12/23
# Resampling test for RNA results

args <- commandArgs(trailingOnly = TRUE)

# load tidyverse
library(tidyverse)

# specify number of permutations
P <- 20000

# read in input file
# store all files in a list for Roller data
fileName <- args[1]
label <- strsplit(tail(strsplit(fileName, '/')[[1]], 1), '_ENSEMBLorthologs.txt')[[1]][1]

# store job number
# jobNum <- args[2]

# GET UNIQUE EXON LENGTHS IN EACH SPECIES
# from a file make in Fig6_RNA.R
exon.lengths <- read_delim('exon.lengths.txt')

### Store individual species pair tables in one giant table

# store file names in a list
filePathToData <- '/home/ak2267/project/Roller/RNA/221121_RNA/ENSEMBL_orthologs/'
#filePathToData <- '~/Desktop/CGI/RNA/ENSEMBL_orthologs/'

# read in table
data <- read_delim(paste0(filePathToData, fileName), delim = '\t')

# separate file name into columns and use left_join to add exon lengths from above
table.int <- left_join(data, exon.lengths %>% dplyr::select(length.bp, ENSG_A), by = 'ENSG_A')
colnames(table.int)[length(colnames(table.int))] <- 'length.bp.A'
table <- left_join(table.int, exon.lengths %>% dplyr::select(length.bp, ENSG_B), by = 'ENSG_B')
colnames(table)[length(colnames(table))] <- 'length.bp.B'

# drop rows with NA values in any of the count columns
# there are NAs when a 1:1 ortholog isn't actually annotated in the GTF file for that species
# affects <1% of lines so I'm proceeding
table <- table %>%
  drop_na(Counts_A1, Counts_A2, Counts_A3,
          Counts_B1, Counts_B2, Counts_B3)

# replace Counts = 0 with 1
table[table$Counts_A1==0,]$Counts_A1 <- 0.01
table[table$Counts_A2==0,]$Counts_A2 <- 0.01
table[table$Counts_A3==0,]$Counts_A3 <- 0.01
table[table$Counts_B1==0,]$Counts_B1 <- 0.01
table[table$Counts_B2==0,]$Counts_B2 <- 0.01
table[table$Counts_B3==0,]$Counts_B3 <- 0.01

# get RPK
table$RPK_A1 <- table$Counts_A1 / table$length.bp.A
table$RPK_A2 <- table$Counts_A2 / table$length.bp.A
table$RPK_A3 <- table$Counts_A3 / table$length.bp.A
table$RPK_B1 <- table$Counts_B1 / table$length.bp.B
table$RPK_B2 <- table$Counts_B2 / table$length.bp.B
table$RPK_B3 <- table$Counts_B3 / table$length.bp.B

# get sum of RPKs for each group
sumOfRPKs <- table %>%
  summarize(RPKsum_A1 = sum(RPK_A1),
            RPKsum_A2 = sum(RPK_A2),
            RPKsum_A3 = sum(RPK_A3),
            RPKsum_B1 = sum(RPK_B1),
            RPKsum_B2 = sum(RPK_B2),
            RPKsum_B3 = sum(RPK_B3))

# calculate TPMs
table$TPM_A1 <- table$RPK_A1 / sumOfRPKs$RPKsum_A1 * (10^6)
table$TPM_A2 <- table$RPK_A2 / sumOfRPKs$RPKsum_A2 * (10^6)
table$TPM_A3 <- table$RPK_A3 / sumOfRPKs$RPKsum_A3 * (10^6)
table$TPM_B1 <- table$RPK_B1 / sumOfRPKs$RPKsum_B1 * (10^6)
table$TPM_B2 <- table$RPK_B2 / sumOfRPKs$RPKsum_B2 * (10^6)
table$TPM_B3 <- table$RPK_B3 / sumOfRPKs$RPKsum_B3 * (10^6)

table$TPM.A <- (table$TPM_A1 + table$TPM_A2 + table$TPM_A3) / 3
table$TPM.B <- (table$TPM_B1 + table$TPM_B2 + table$TPM_B3) / 3

table$TPM.ratio <- table$TPM.A / table$TPM.B

# print median OBSERVED TPM.ratio for A and B to observed file
median_obs_A <- table %>% filter(CGI_summary == 'A') %>% select(TPM.ratio) %>% summarize(median = median(TPM.ratio))
median_obs_B <- table %>% filter(CGI_summary == 'B') %>% select(TPM.ratio) %>% summarize(median = median(TPM.ratio))
write.table(c(median_obs_A, median_obs_B), file = paste0(label, '_observedMedians.txt'), quote = F, row.names = F, col.names = F)

###################################################
## EXPANDED NA.medians PROCEDURE WITH RESAMPLING ##
###################################################

# resample ratio across all TPM level bins
# since ratio can be skewed depending on expression level

# calculate geometric mean of TPMs between speciesA and speciesB
table$TPM_Avg <- sqrt( table$TPM.A * table$TPM.B)

# store resampling parameters
N = P
binSizeToCutoff1 <- 1
cutoff1 <- 10
binSizeToCutoff2 <- 10
cutoff2 <- 500
binSizeAfterCutoff2 <- 200

# take A & B in subTable and NA in backgroundTable
backgroundTable <- table %>% filter(CGI_summary %in% c(NA))
  
# initiate resampleTable (will add a row with values every resample)
resampleList_A <- c()
resampleList_B <- c()

for (i in 1:N) {
    
  # FIRST DO A-ONLY
  sampleList <- c()
  subTable <- table %>% filter(CGI_summary == 'A')
    
  maxTPM <- max(subTable$TPM_Avg)
  binNumberToCutoff1 <- ceiling( cutoff1 / binSizeToCutoff1 )
  binNumberToCutoff2 <- ceiling( (cutoff2 - cutoff1) / binSizeToCutoff2 )
  binNumberAfterCutoff2 <- ceiling( (maxTPM - cutoff2) / binSizeAfterCutoff2 )
    
  for (j in 1:binNumberToCutoff1) {
    #print(j)
    binLeft <- binSizeToCutoff1 * (j-1)
    binRight <- binSizeToCutoff1 * j
      
    n_inBin <- nrow(subTable[subTable$TPM_Avg >= binLeft & subTable$TPM_Avg < binRight,])
    
    if (n_inBin > 0) {
      backgroundTableInBin <- backgroundTable[backgroundTable$TPM_Avg >= binLeft & backgroundTable$TPM_Avg < binRight,]
      
      while (nrow(backgroundTableInBin) == 0) {
        binLeft <- binLeft - binSizeToCutoff1
        binRight <- binRight + binSizeToCutoff1
        backgroundTableInBin <- backgroundTable[backgroundTable$TPM_Avg >= binLeft & backgroundTable$TPM_Avg < binRight,]
      }
      sampleFrom <- as.vector(backgroundTableInBin$TPM.ratio)
      sampleList <- c(sampleList, sampleFrom[sample(1:nrow(backgroundTableInBin), n_inBin, replace = T)])
    }
  }
    
  for (j in 1:binNumberToCutoff2) {
    #print(j)
    binLeft <- cutoff1 + binSizeToCutoff2 * (j-1)
    binRight <- cutoff1 + binSizeToCutoff2 * j
      
    n_inBin <- nrow(subTable[subTable$TPM_Avg >= binLeft & subTable$TPM_Avg < binRight,])
    
    if (n_inBin > 0) {
      backgroundTableInBin <- backgroundTable[backgroundTable$TPM_Avg >= binLeft & backgroundTable$TPM_Avg < binRight,]
      
      while (nrow(backgroundTableInBin) == 0) {
        binLeft <- binLeft - binSizeToCutoff2
        binRight <- binRight + binSizeToCutoff2
        backgroundTableInBin <- backgroundTable[backgroundTable$TPM_Avg >= binLeft & backgroundTable$TPM_Avg < binRight,]
      }
      sampleFrom <- as.vector(backgroundTableInBin$TPM.ratio)
      sampleList <- c(sampleList, sampleFrom[sample(1:nrow(backgroundTableInBin), n_inBin, replace = T)])
    }
  }
    
  for (j in 1:binNumberAfterCutoff2) {
    #print(j)
    binLeft <- cutoff2 + binSizeAfterCutoff2 * (j-1)
    binRight <- cutoff2 + binSizeAfterCutoff2 * j
    
    n_inBin <- nrow(subTable[subTable$TPM_Avg >= binLeft & subTable$TPM_Avg < binRight,])
    
    if (n_inBin > 0) {
      backgroundTableInBin <- backgroundTable[backgroundTable$TPM_Avg >= binLeft & backgroundTable$TPM_Avg < binRight,]
      
      while (nrow(backgroundTableInBin) == 0) {
        binLeft <- binLeft - binSizeAfterCutoff2
        binRight <- binRight + binSizeAfterCutoff2
        backgroundTableInBin <- backgroundTable[backgroundTable$TPM_Avg >= binLeft & backgroundTable$TPM_Avg < binRight,]
      }
      sampleFrom <- as.vector(backgroundTableInBin$TPM.ratio)
      sampleList <- c(sampleList, sampleFrom[sample(1:nrow(backgroundTableInBin), n_inBin, replace = T)])
    }
  }  
  
  # median from sample
  median <- median(sampleList, na.rm = T)
    
  # add median to resampleList
  resampleList_A <- c(resampleList_A, median)
    
  # THEN DO B-ONLY
  sampleList <- c()
  subTable <- table %>% filter(CGI_summary == 'B')
      
  maxTPM <- max(subTable$TPM_Avg)
  binNumberToCutoff1 <- ceiling( cutoff1 / binSizeToCutoff1 )
  binNumberToCutoff2 <- ceiling( (cutoff2 - cutoff1) / binSizeToCutoff2 )
  binNumberAfterCutoff2 <- ceiling( (maxTPM - cutoff2) / binSizeAfterCutoff2 )
      
  for (j in 1:binNumberToCutoff1) {
    #print(j)
    binLeft <- binSizeToCutoff1 * (j-1)
    binRight <- binSizeToCutoff1 * j
        
    n_inBin <- nrow(subTable[subTable$TPM_Avg >= binLeft & subTable$TPM_Avg < binRight,])
    
    if (n_inBin > 0) {
      backgroundTableInBin <- backgroundTable[backgroundTable$TPM_Avg >= binLeft & backgroundTable$TPM_Avg < binRight,]
      
      while (nrow(backgroundTableInBin) == 0) {
        binLeft <- binLeft - binSizeToCutoff1
        binRight <- binRight + binSizeToCutoff1
        backgroundTableInBin <- backgroundTable[backgroundTable$TPM_Avg >= binLeft & backgroundTable$TPM_Avg < binRight,]
      }
      sampleFrom <- as.vector(backgroundTableInBin$TPM.ratio)
      sampleList <- c(sampleList, sampleFrom[sample(1:nrow(backgroundTableInBin), n_inBin, replace = T)])
    }
  }
 
  for (j in 1:binNumberToCutoff2) {
    #print(j)
    binLeft <- cutoff1 + binSizeToCutoff2 * (j-1)
    binRight <- cutoff1 + binSizeToCutoff2 * j
        
    n_inBin <- nrow(subTable[subTable$TPM_Avg >= binLeft & subTable$TPM_Avg < binRight,])
    
    if (n_inBin > 0) {
      backgroundTableInBin <- backgroundTable[backgroundTable$TPM_Avg >= binLeft & backgroundTable$TPM_Avg < binRight,]
      
      while (nrow(backgroundTableInBin) == 0) {
        binLeft <- binLeft - binSizeToCutoff2
        binRight <- binRight + binSizeToCutoff2
        backgroundTableInBin <- backgroundTable[backgroundTable$TPM_Avg >= binLeft & backgroundTable$TPM_Avg < binRight,]
      }
      sampleFrom <- as.vector(backgroundTableInBin$TPM.ratio)
      sampleList <- c(sampleList, sampleFrom[sample(1:nrow(backgroundTableInBin), n_inBin, replace = T)])
    }
  }
     
  for (j in 1:binNumberAfterCutoff2) {
    #print(j)
    binLeft <- cutoff2 + binSizeAfterCutoff2 * (j-1)
    binRight <- cutoff2 + binSizeAfterCutoff2 * j
        
    n_inBin <- nrow(subTable[subTable$TPM_Avg >= binLeft & subTable$TPM_Avg < binRight,])
    
    if (n_inBin > 0) {
      backgroundTableInBin <- backgroundTable[backgroundTable$TPM_Avg >= binLeft & backgroundTable$TPM_Avg < binRight,]
      
      while (nrow(backgroundTableInBin) == 0) {
        binLeft <- binLeft - binSizeAfterCutoff2
        binRight <- binRight + binSizeAfterCutoff2
        backgroundTableInBin <- backgroundTable[backgroundTable$TPM_Avg >= binLeft & backgroundTable$TPM_Avg < binRight,]
      }
      sampleFrom <- as.vector(backgroundTableInBin$TPM.ratio)
      sampleList <- c(sampleList, sampleFrom[sample(1:nrow(backgroundTableInBin), n_inBin, replace = T)])
    }
    
  }  
  # median from sample
  median <- median(sampleList, na.rm = T)
      
  # add median to resampleList
  resampleList_B <- c(resampleList_B, median)

}

# write results (resampleList_A and resampleList_B) to output files for use in summarizing downstream
write.table(resampleList_A, paste0(label, '_A', '_resamplingMedians.txt'), sep = '\n', quote = F, row.names = F, col.names = F)
write.table(resampleList_B, paste0(label, '_B', '_resamplingMedians.txt'), sep = '\n', quote = F, row.names = F, col.names = F)

