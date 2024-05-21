# 10/24/22
# Purpose: generate histogram for example permutation test in Fig S22 (originally S17)

require(ggplot2)
require(cowplot)
theme_set(theme_cowplot())

setwd('/Users/acadiak/Desktop/CGI/Figures/Mar_2023')

### modify permutation.test to output counts and distribution for example - rheMac10 vs mm39 brain H3K27ac
label <- 'rheMac10_mm39_brain_H3K4me3'

# define order of CGI status
CGI_order <- c('A','B','both')
# define order of Peak status
Peak_order <- c('A','B','both','0')

# read in data
filePathRoller <- '/Users/acadiak/Desktop/CGI/speciesPairs/CGIcentric/Roller_summaryFiles/'
file <- paste(filePathRoller,label,'.txt',sep='')
table <- read.table(file,header=T)

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
subsetTable_3x3 <- table[table$Peak_A==1|table$Peak_B==1,]

# reduce to only necessary columns
subsetTable_3x3 <- subsetTable_3x3[,c(39,40)]
inputTable <- subsetTable_3x3

# get sample info from label
speciesA <- strsplit(label,'_')[[1]][1]
speciesB <- strsplit(label,'_')[[1]][2]
speciesPair <- paste(speciesA,'_',speciesB,sep='')
tissue <- strsplit(label,'_')[[1]][3]
mark <- strsplit(label,'_')[[1]][4]
infoList <- c(speciesPair, speciesA, speciesB, tissue, mark)

# collect all possible CGI x Peak combos
comboList <- c()
for (CGI in CGI_order) {
  for (Peak in Peak_order[1:3]) {
    comboList <- c(comboList, paste(CGI,Peak,sep='_'))
  }
}

# set observed values
observedValues <- c()
for (CGI in CGI_order) {
  for (Peak in Peak_order[1:3]) {
    observedValues <- c(observedValues, nrow(inputTable[inputTable$CGI_summary==CGI&
                                                          inputTable$Peak_summary==Peak,]))
  }
}

# calculate expected values
dimensions = 3
observedTable <- matrix(observedValues, nrow = dimensions)
expectedTable <- observedTable
for (i in c(1,2,3)) {
  for (j in seq(1,dimensions,1)) {
    expectedTable[j,i] <- (sum(observedTable[j,]) * sum(observedTable[,i])) / sum(observedTable)
  }
}
expectedValues <- as.vector(expectedTable)


#########################
# RUN PERMUTATION TEST #
#########################

# for making histogram in Fig S22 (explanation of permutation test)

# prepare table to collect numbers in each CGI_Peak category
P <- 20000
set.seed(7)

countTable <- matrix(nrow = P,
                     ncol = length(comboList))
colnames(countTable) <- comboList

# perform permutation
for (i in 1:P) {
  # permute Peak category
  inputTable$permuted <- sample(inputTable$Peak_summary)

  # count number of sites in each CGI_Peak category
  for (combo in comboList) {
    CGI <- strsplit(combo, '_')[[1]][1]
    Peak <- strsplit(combo, '_')[[1]][2]
    countTable[i, combo] <- nrow(inputTable[inputTable$CGI_summary==CGI&
                                              inputTable$permuted==Peak,])
  }
}

# view results for values in tables in Fig S22
observedTable
round(expectedTable, 0)
round(log2(observedTable / expectedTable), 2)

# observedTable
# 163   21   39
# 15  143   31
# 11    5   29
# expectedTable
# 92   82   48
# 78   70   41
# 19   17   10
# log2(observedTable / expectedTable)
# 0.82 -1.97 -0.31
# -2.38  1.03 -0.40
# -0.76 -1.73  1.57

# plot histogram for Fig S22
countTable.df <- data.frame(countTable)

pdf(file = 'FigS17_exampleHistogram.pdf', height = 3, width = 4)
hist(countTable.df$A_A, 
     breaks = 20,
     xlim = c(60, 180),
     ylim = c(0, 3000),
     las = 1,
     xlab = '',
     main = '')
abline(v = expectedValues[1], col = 'black', lwd = 2)
abline(v = observedValues[1], col = 'red', lwd = 2)
dev.off()

### see main Fig3 script for heatmap generation

