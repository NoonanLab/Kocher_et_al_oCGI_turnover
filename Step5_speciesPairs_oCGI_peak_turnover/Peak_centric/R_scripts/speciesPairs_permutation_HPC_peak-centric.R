# 1/9/22
# Purpose: implement permutation test for speciesPair analysis in parallel jobs on Ruddle
# write Fold Difference and p values to a text file for use downstream
# This script is for PEAK-CENTRIC analysis, as opposed to CGI-centric

args <- commandArgs(trailingOnly = TRUE)

# set output directory
setwd('/gpfs/gibbs/pi/noonan/ak2267/speciesPairs/peakCentric/permutation')

# specify number of permutations
P <- 20000

# read in input file
# store all files in a list for Roller data
filePath <- args[1]
label <- strsplit(tail(strsplit(filePath, '/')[[1]], 1), '.txt')[[1]][1]

# define order of Peak status
Peak_order <- c('A','B','both')
# define order of CGI status
CGI_order <- c('A','B','both','0')

# get sample info from label
speciesA <- strsplit(label,'_')[[1]][1]
speciesB <- strsplit(label,'_')[[1]][2]
speciesPair <- paste(speciesA,'_',speciesB,sep='')
tissue <- strsplit(label,'_')[[1]][3]
mark <- strsplit(label,'_')[[1]][4]
tissueTranslator <- data.frame(shorthand=c('brain', 'limb'),
                                  longhand=c('devBrain', 'devLimb'))
markTranslator <- data.frame(shorthand=c('ac', 'me2'),
                             longhand=c('H3K27ac', 'H3K4me2'))
if (speciesA=='hg19' || speciesA=='rheMac2') {
  timePoint <- strsplit(label, '_')[[1]][5]
  tissue <- tissueTranslator[tissueTranslator$shorthand==tissue,]$longhand
  mark <- markTranslator[markTranslator$shorthand==mark,]$longhand
}
if (speciesA!='hg19' && speciesA!='rheMac2') {
  timePoint <- NA
}

# read in data
table <- read.table(filePath, header=T)

# add column describing species-specificity of CGI
table$CGI_summary <- '0'
table[table$CGI_A>0&table$CGI_B>0,]$CGI_summary <- 'both'
table[table$CGI_A>0&table$CGI_B==0,]$CGI_summary <- 'A'
table[table$CGI_A==0&table$CGI_B>0,]$CGI_summary <- 'B'

# add column describing species-specificity of peak
table$Peak_summary <- '0'
table[table$Peak_A>0&table$Peak_B>0,]$Peak_summary <- 'both'
table[table$Peak_A>0&table$Peak_B==0,]$Peak_summary <- 'A'
table[table$Peak_A==0&table$Peak_B>0,]$Peak_summary <- 'B'

# work with either 3 x 3 table (none without activity) or 4 x 3 table (full dataset)
subsetTable_3x3 <- table[table$CGI_A==1|table$CGI_B==1,]
subsetTable_4x3 <- table

# reduce to only necessary columns
subsetTable_3x3 <- subsetTable_3x3[,c('Peak_summary','CGI_summary')]
subsetTable_4x3 <- subsetTable_4x3[,c('Peak_summary','CGI_summary')]

# get sample info from label
infoList <- c(speciesPair, speciesA, speciesB, tissue, mark, timePoint)

###### run permutation test with 3x3 table
dimensions <- 3
inputTable <- subsetTable_3x3

# collect all possible CGI x Peak combos
comboList <- c()
for (Peak in Peak_order) {
  for (CGI in CGI_order[1:dimensions]) {
    comboList <- c(comboList, paste(Peak,CGI,sep='_'))
    }
}
  
# set observed values
observedValues <- c()
for (Peak in Peak_order) {
  for (CGI in CGI_order[1:dimensions]) {
    observedValues <- c(observedValues, nrow(inputTable[inputTable$Peak_summary==Peak&
                                                          inputTable$CGI_summary==CGI,]))
  }
}

# calculate expectedValues (row sum * col sum / total sum)
# in table where columns are Peaks (A, B, both) and rows are CGI status (A, B, both, <neither>)
observedTable <- matrix(observedValues, nrow = dimensions)
expectedTable <- observedTable
for (i in c(1,2,3)) {
  for (j in seq(1,dimensions,1)) {
    expectedTable[j,i] <- (sum(observedTable[j,]) * sum(observedTable[,i])) / sum(observedTable)
  }
}
expectedValues <- as.vector(expectedTable)

# prepare table to collect numbers in each CGI_Peak category
countTable <- matrix(nrow = P,
                     ncol = length(comboList))
colnames(countTable) <- comboList

# perform permutation
for (i in 1:P) {
  # permute Peak category
  inputTable$permuted <- NA
  inputTable$permuted <- sample(inputTable$CGI_summary)
  
  # count number of sites in each Peak_CGI category
  tableSummary <- data.frame(table(inputTable[,-2]))
  tableSummary$name <- paste(tableSummary$Peak_summary, tableSummary$permuted, sep = '_')
  tableSummary$index <- match(comboList, tableSummary$name)
  countTable[i, ] <- tableSummary[order(tableSummary$index),]$Freq
}

# calculate FD and p value
FD_list <- c()
p_list <- c()
for (i in 1:length(comboList)) {
  FD_list <- c(FD_list, observedValues[i] / expectedValues[i])
  p_list <- c(p_list, max(min(length(countTable[,i][countTable[,i] < observedValues[i]]),
                          length(countTable[,i][countTable[,i] > observedValues[i]])), 1) / P)
}

# return label, FD, and p value
write.table(t(data.frame(c(infoList, FD_list, p_list))), 
            file = paste(label, '_', dimensions, 'x3.txt', sep = ''), 
            sep = '\t', row.names = F, col.names = F, quote = F)

###### run permutation test with 4x3 table
dimensions <- 4
inputTable <- subsetTable_4x3

# collect all possible CGI x Peak combos
comboList <- c()
for (Peak in Peak_order) {
  for (CGI in CGI_order[1:dimensions]) {
    comboList <- c(comboList, paste(Peak,CGI,sep='_'))
  }
}

# set observed values
observedValues <- c()
for (Peak in Peak_order) {
  for (CGI in CGI_order[1:dimensions]) {
    observedValues <- c(observedValues, nrow(inputTable[inputTable$Peak_summary==Peak&
                                                          inputTable$CGI_summary==CGI,]))
  }
}

# calculate expectedValues (row sum * col sum / total sum)
# in table where columns are Peaks (A, B, both) and rows are CGI status (A, B, both, <neither>)
observedTable <- matrix(observedValues, nrow = dimensions)
expectedTable <- observedTable
for (i in c(1,2,3)) {
  for (j in seq(1,dimensions,1)) {
    expectedTable[j,i] <- (sum(expectedTable[j,]) * sum(expectedTable[,i])) / sum(expectedTable)
  }
}
expectedValues <- as.vector(expectedTable)

# prepare table to collect numbers in each CGI_Peak category
countTable <- matrix(nrow = P,
                     ncol = length(comboList))
colnames(countTable) <- comboList

# perform permutation
for (i in 1:P) {
  # permute Peak category
  inputTable$permuted <- NA
  inputTable$permuted <- sample(inputTable$CGI_summary)
  
  # count number of sites in each Peak_CGI category
  tableSummary <- data.frame(table(inputTable[,-2]))
  tableSummary$name <- paste(tableSummary$Peak_summary, tableSummary$permuted, sep = '_')
  tableSummary$index <- match(comboList, tableSummary$name)
  countTable[i, ] <- tableSummary[order(tableSummary$index),]$Freq
}

# calculate FD and p value
FD_list <- c()
p_list <- c()
for (i in 1:length(comboList)) {
  FD_list <- c(FD_list, observedValues[i] / expectedValues[i])
  p_list <- c(p_list, max(min(length(countTable[,i][countTable[,i] < observedValues[i]]),
                              length(countTable[,i][countTable[,i] > observedValues[i]])), 1) / P)
}

# return label, FD, and p value
write.table(t(data.frame(c(infoList, FD_list, p_list))), 
            file = paste(label, '_', dimensions, 'x3.txt', sep = ''), 
            sep = '\t', row.names = F, col.names = F, quote = F)
