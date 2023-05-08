# 2/15/23
# Purpose: run resampling test for enrichment of HGEs with specific oCGI patterns

args <- commandArgs(trailingOnly = TRUE)

# load tidyverse
library(tidyverse)

# specify number of permutations
P <- 20000

# read in input file
# store all files in a list for Roller data
fileInfo <- args[1]
fileName <- paste0('HGEtable_', fileInfo, '.txt')

filePathToData <- '/home/ak2267/project/EnhancerClasses/HGE/HGEtable/'

# read in table
data <- read_delim(paste0(filePathToData, fileName), delim = '\t')

# filter based on being an enhancer and lifting to rhesus and mouse
table <- data %>%
  filter(Enhancer_AK == 1,
         LiftsToRhesus==1,
         LiftsToMouse==1)

# add average signals
table$avgSignal = (table$avgSignal_rep1 + table$avgSignal_rep2) / 2
table$totalSignal = (table$totalSignal_rep1 + table$totalSignal_rep2) / 2

# add column with CGI pattern, i.e. H, HR, HRM, RM, etc
table$CGI_pattern <- 'None'
table[table$CGI_human == 1 & table$CGI_rhesus == 1 & table$CGI_mouse == 1,]$CGI_pattern <- 'HRM'
table[table$CGI_human == 1 & table$CGI_rhesus == 1 & table$CGI_mouse == 0,]$CGI_pattern <- 'HR'
table[table$CGI_human == 1 & table$CGI_rhesus == 0 & table$CGI_mouse == 1,]$CGI_pattern <- 'HM'
table[table$CGI_human == 0 & table$CGI_rhesus == 1 & table$CGI_mouse == 1,]$CGI_pattern <- 'RM'
table[table$CGI_human == 1 & table$CGI_rhesus == 0 & table$CGI_mouse == 0,]$CGI_pattern <- 'H'
table[table$CGI_human == 0 & table$CGI_rhesus == 1 & table$CGI_mouse == 0,]$CGI_pattern <- 'R'
table[table$CGI_human == 0 & table$CGI_rhesus == 0 & table$CGI_mouse == 1,]$CGI_pattern <- 'M'

###################
# Resampling test #
###################

# set binSize. with avgSignal, would use binSize = 0.05
# but here we'll use totalSignal
# two steps:
binSizeToCutoff <- 10
cutoff <- 2000
binSizeAfterCutoff <- 100

# get foreground (HGEs) and background (non-HGEs)
gainTable <- table %>% filter(HGE == 1)
nonGainTable <- table %>% filter(HGE == 0)
  
# get max signal in foreground
maxSignal <- max(gainTable$totalSignal)

# calculate number of bins
binNumberToCutoff <- ceiling( cutoff / binSizeToCutoff )
binNumberAfterCutoff <- ceiling( (maxSignal - cutoff) / binSizeAfterCutoff )

resultsTable <- tibble(
  'None' = numeric(),
  'H' = numeric(),
  'R' = numeric(),
  'M' = numeric(),
  'HR' = numeric(),
  'HM' = numeric(),
  'RM' = numeric(),
  'HRM' = numeric(),
)

for (i in 1:P) {
  
  # initialize data frame to store sample info
  resampleTable <- data.frame()
  
  # sample from each bin
  for (j in 1:binNumberToCutoff) {
      
    binLeft <- binSizeToCutoff * (j-1)
    binRight <- binSizeToCutoff * j
      
    numberGains <- nrow(gainTable[gainTable$totalSignal>binLeft & gainTable$totalSignal<=binRight,])
    
    if (numberGains > 0) {
      nonGainsInTheBin <- nonGainTable[nonGainTable$totalSignal>binLeft & nonGainTable$totalSignal<=binRight,]

      while (nrow(nonGainsInTheBin) == 0) {
        binLeft <- binLeft - binSizeToCutoff
        binRight <- binRight + binSizeToCutoff
        nonGainsInTheBin <- nonGainTable[nonGainTable$totalSignal>binLeft & nonGainTable$totalSignal<=binRight,]
      }
      sample <- nonGainsInTheBin[sample(1:nrow(nonGainsInTheBin), numberGains, replace = T),]
      
      # append the sample to resampleTable
      resampleTable <- rbind(resampleTable, sample)
    }
      
  }
  
  for (k in 1:binNumberAfterCutoff) {
    binLeft <- cutoff + binNumberAfterCutoff * (k-1)
    binRight <- cutoff + binNumberAfterCutoff * k
    
    numberGains <- nrow(gainTable[gainTable$totalSignal>binLeft & gainTable$totalSignal<=binRight,])
    
    if (numberGains > 0) {
      nonGainsInTheBin <- nonGainTable[nonGainTable$totalSignal>binLeft & nonGainTable$totalSignal<=binRight,]
      
      while (nrow(nonGainsInTheBin) == 0) {
        binLeft <- binLeft - binSizeAfterCutoff
        binRight <- binRight + binSizeAfterCutoff
        nonGainsInTheBin <- nonGainTable[nonGainTable$totalSignal>binLeft & nonGainTable$totalSignal<=binRight,]
      }
      
      sample <- nonGainsInTheBin[sample(1:nrow(nonGainsInTheBin), numberGains, replace = T),]
      
      # append the sample to resampleTable
      resampleTable <- rbind(resampleTable, sample)
    }
    
  }
  
  # do counting for this P and append to resultsTable
  counts <- resampleTable %>% group_by(CGI_pattern) %>% count() %>% pivot_wider(names_from = CGI_pattern, values_from = n)
  resultsTable <- full_join(resultsTable, counts)
  
}


# write observed values within HGEs to an output table
observed_HGEs <- gainTable %>% group_by(CGI_pattern) %>% count() %>% pivot_wider(names_from = CGI_pattern, values_from = n)
write.table(observed_HGEs, paste0(fileInfo, '_observed_HGEs_CGI_patterns.txt'), sep = '\t', quote = F, row.names = F, col.names = T)

# write results of resampling 
write.table(resultsTable, paste0(fileInfo, '_resampled_HGEs_CGI_patterns.txt'), sep = '\t', quote = F, row.names = F, col.names = T)

