# 11/22/22
# Purpose: analyze + make plots for Roller RNA data

require(tidyverse)
require(RColorBrewer)
require(patchwork)

# change this later to a figure folder
setwd('/Users/acadiak/Desktop/CGI/Figures/Mar_2023')

# GET UNIQUE EXON LENGTHS IN EACH SPECIES

# package to get exon lengths
# https://www.biostars.org/p/83901/
BiocManager::install("GenomicFeatures")
library('GenomicFeatures')

GTFpath <- '/Users/acadiak/Desktop/CGI/RNA/GTF/'

speciesList <- c('rheMac10', 'mm39', 'rn7', 'susScr11', 'felCat9', 'equCab3')
exon.lengths <- tibble()
for (species in speciesList) {
  txdb <- makeTxDbFromGFF(paste0(GTFpath, species, '.ENSEMBLinUCSC.gtf'), format="gtf")
  exons.list.per.gene <- exonsBy(txdb, by="gene")
  exonic.gene.sizes <- as.data.frame(sum(width(GenomicRanges::reduce(exons.list.per.gene))))
  exon.lengths <- bind_rows(exon.lengths, exonic.gene.sizes)
}
colnames(exon.lengths) <- 'length.bp'
exon.lengths$ENSG_A <- rownames(exon.lengths)
exon.lengths$ENSG_B <- exon.lengths$ENSG_A

# write exon.lengths to a text file and upload to the cluster for use in resampling jobs
write_tsv(exon.lengths, file = '/Users/acadiak/Desktop/CGI/RNA/exon.lengths.txt', quote = 'none')
exon.lengths <- read_tsv('/Users/acadiak/Desktop/CGI/RNA/exon.lengths.txt') %>% tibble()

######## Read in data

# store file names in a list
filePathToData <- '/Users/acadiak/Desktop/CGI/RNA/ENSEMBL_orthologs'
files <- dir(path = filePathToData, pattern = "*_ENSEMBLorthologs.txt")

# read in as nested dataframe
data <- tibble(filename = files) %>%
  mutate(file_contents = map(filename,
                             function(x) read_tsv(file.path(filePathToData, x)))
  )

# unnest to generate full table, but with file name as a column
unnested <- unnest(data, cols = c(file_contents)) %>%
  separate(filename,
           c('speciesA', 'speciesB', 'tissue', 'mark'),
           fill = 'right')
unnested$speciesPair <- paste(unnested$speciesA, unnested$speciesB, sep = '_')
unnested$group <- paste(unnested$speciesPair, unnested$tissue, unnested$mark, sep = '_')

# separate file name into columns and use left_join to add exon lengths from above
table.int <- left_join(unnested, exon.lengths %>% dplyr::select(length.bp, ENSG_A), by = 'ENSG_A')
colnames(table.int)[length(colnames(table.int))] <- 'length.bp.A'
table <- left_join(table.int, exon.lengths %>% dplyr::select(length.bp, ENSG_B), by = 'ENSG_B')
colnames(table)[length(colnames(table))] <- 'length.bp.B'

# drop rows with NA values in any of the count columns
table <- table %>%
  drop_na(Counts_A1, Counts_A2, Counts_A3,
          Counts_B1, Counts_B2, Counts_B3)

# replace Counts = 0 with 0.1
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
  group_by(group) %>% 
  summarize(RPKsum_A1 = sum(RPK_A1),
            RPKsum_A2 = sum(RPK_A2),
            RPKsum_A3 = sum(RPK_A3),
            RPKsum_B1 = sum(RPK_B1),
            RPKsum_B2 = sum(RPK_B2),
            RPKsum_B3 = sum(RPK_B3))
table <- left_join(table, sumOfRPKs, by = 'group')

# calculate TPMs
table$TPM_A1 <- table$RPK_A1 / table$RPKsum_A1 * (10^6)
table$TPM_A2 <- table$RPK_A2 / table$RPKsum_A2 * (10^6)
table$TPM_A3 <- table$RPK_A3 / table$RPKsum_A3 * (10^6)
table$TPM_B1 <- table$RPK_B1 / table$RPKsum_B1 * (10^6)
table$TPM_B2 <- table$RPK_B2 / table$RPKsum_B2 * (10^6)
table$TPM_B3 <- table$RPK_B3 / table$RPKsum_B3 * (10^6)

table$TPM.A <- (table$TPM_A1 + table$TPM_A2 + table$TPM_A3) / 3
table$TPM.B <- (table$TPM_B1 + table$TPM_B2 + table$TPM_B3) / 3

table$TPM.ratio <- table$TPM.A / table$TPM.B

######## Run resampling test - see 221115_RollerRNA.sh, which uses RNA_resamplingTest.R

######## Read results of resampling test

RNAtable <- read.table('~/Desktop/CGI/RNA/RNA_resamplingSummary.txt', header = F) %>% tibble()
colnames(RNAtable) <- c('speciesPair', 'tissue', 'mark', 'CGI', 'observed', 'expected', 'p')
RNAtable$group <- paste(RNAtable$speciesPair, RNAtable$tissue, RNAtable$mark, sep = '_')

# left_join these expected values into the full table for normalization
table <- left_join(table, RNAtable %>% dplyr::select(group, expected), by = 'group')

# normalize TPM.ratio by expected ratio from resampling test
table$TPM.ratio.norm.resampling <- table$TPM.ratio / table$expected

######### PLOTS

# set up p values with x and y values plotting values and stars
RNAtable$padj <- p.adjust(RNAtable$p, method = 'BH')
RNAtable$x <- 1
RNAtable[RNAtable$CGI=='B',]$x <- 2
RNAtable$y <- 1.6
RNAtable$stars <- ''
RNAtable[RNAtable$padj<0.05,]$stars <- '*'

# change species names
RNAtable <- RNAtable %>%
  separate(speciesPair, into = c('speciesA','speciesB')) %>%
  mutate(speciesA = factor(speciesA,
                           levels = c('rheMac10', 'mm39', 'rn7', 'susScr11', 'felCat9'),
                           labels = c('Rhesus', 'Mouse', 'Rat', 'Pig', 'Cat')),
         speciesB = factor(speciesB,
                           levels = c('mm39', 'rn7', 'susScr11', 'felCat9', 'equCab3'),
                           labels = c('Mouse', 'Rat', 'Pig', 'Cat', 'Horse')),
         tissue = factor(tissue,
                         levels = c('brain', 'liver', 'muscle', 'testis'),
                         labels = c('Brain', 'Liver', 'Muscle', 'Testis')))
RNAtable$speciesPair <- paste(RNAtable$speciesA, RNAtable$speciesB, sep = '\n')

# summarized data
summaryTable <- table %>%
  filter(CGI_summary %in% c('A', 'B')) %>%
  mutate(speciesA = factor(speciesA,
                           levels = c('rheMac10', 'mm39', 'rn7', 'susScr11', 'felCat9'),
                           labels = c('Rhesus', 'Mouse', 'Rat', 'Pig', 'Cat')),
         speciesB = factor(speciesB,
                           levels = c('mm39', 'rn7', 'susScr11', 'felCat9', 'equCab3'),
                           labels = c('Mouse', 'Rat', 'Pig', 'Cat', 'Horse'))) %>%
  group_by(speciesA, speciesB, tissue, mark, CGI_summary) %>% 
  summarise(Q1=quantile(log2(TPM.ratio.norm.resampling), probs = 0.25),
            median=median(log2(TPM.ratio.norm.resampling)), 
            Q3=quantile(log2(TPM.ratio.norm.resampling), probs = 0.75),
            n.obs = n())
summaryTable$color <- 'black'
summaryTable[summaryTable$CGI_summary=='A',]$color <- 'royalblue4'
summaryTable[summaryTable$CGI_summary=='B',]$color <- 'lightskyblue'
summaryTable$speciesPair <- paste(summaryTable$speciesA, summaryTable$speciesB, sep = '\n')


# PLOTS

summaryTable$tissue <- factor(summaryTable$tissue,
                              levels = c('brain', 'liver', 'muscle', 'testis'),
                              labels = c('Brain', 'Liver', 'Muscle', 'Testis'))
summaryTable$mark <- factor(summaryTable$mark,
                            levels = c('H3K27ac', 'H3K4me1', 'H3K4me3'))
summaryTable$color <- factor(summaryTable$color,
                            levels = c('royalblue4', 'lightskyblue'))
summaryTable$speciesPair <- factor(summaryTable$speciesPair,
                                   levels = c('Rhesus\nMouse', 'Rhesus\nRat', 'Rhesus\nPig', 'Rhesus\nCat', 'Rhesus\nHorse',
                                              'Mouse\nRat', 'Mouse\nPig', 'Mouse\nCat', 'Mouse\nHorse',
                                              'Rat\nPig', 'Rat\nCat', 'Rat\nHorse',
                                              'Pig\nCat', 'Pig\nHorse',
                                              'Cat\nHorse'))

function_RNAplots <- function(plotTable, inMark, pointSize) {
  x <- plotTable %>%
    filter(mark==inMark) %>%
    ggplot() +
    geom_hline(yintercept = 0) +
    facet_grid(factor(tissue) ~ factor(speciesPair)) +
    geom_segment(aes(x = CGI_summary, xend = CGI_summary, y = Q1, yend = Q3)) +
    geom_point(aes(x = CGI_summary, y = median, color = color), size = pointSize) +
    scale_color_manual(values = c('royalblue4', 'lightskyblue')) +
    theme_bw() + 
    theme(panel.grid.minor.x = element_blank(),
          panel.grid.major.x = element_blank(),
          panel.grid.minor.y = element_blank(),
          strip.background = element_blank(),
          axis.text.y = element_text(size = 12),
          axis.text.x = element_text(size = 12),
          axis.title = element_text(size = 12),
          strip.text.x = element_text(size = 9),
          strip.text.y = element_text(size = 12),
          panel.spacing.y = unit(24, 'pt'),
          legend.position = 'none') +
    labs(y = 'Log2(TPMA / TPMB)',
         x = 'Species with oCGI and Peak') +
    geom_text(data = RNAtable %>%
                filter(mark==inMark),
              aes(x = x , y = y, label = stars), inherit.aes = F, size = 6)
  return(x)
}

H3K27ac <- function_RNAplots(summaryTable, 'H3K27ac', 2)
H3K4me1 <- function_RNAplots(summaryTable, 'H3K4me1', 2)
H3K4me3 <- function_RNAplots(summaryTable, 'H3K4me3', 2)

# scales = 'free' only for H3K4me3
allTissues <- H3K27ac + H3K4me1 + H3K4me3 + plot_layout(nrow = 3)
ggsave('FigS6.2_RNA_allDatasets.pdf', allTissues, height = 6000, width = 2200, units = 'px')


# make main figure
Fig6 <- summaryTable %>%
  filter(speciesPair %in% c('Rhesus\nRat',
                            'Mouse\nPig',
                            'Rat\nPig',
                            'Pig\nHorse',
                            'Cat\nHorse'),
         tissue == 'Brain') %>%
  ggplot() +
  geom_hline(yintercept = 0) +
  facet_grid(factor(mark) ~ factor(speciesPair), scales = 'free') +
  geom_segment(aes(x = CGI_summary, xend = CGI_summary, y = Q1, yend = Q3)) +
  geom_point(aes(x = CGI_summary, y = median, color = color), size = 2) +
  scale_color_manual(values = c('royalblue4', 'lightskyblue')) +
  theme_bw() + 
  theme(panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        strip.background = element_blank(),
        axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12),
        axis.title = element_text(size = 12),
        strip.text.x = element_text(size = 9),
        strip.text.y = element_text(size = 12),
        legend.position = 'none') +
  labs(y = 'Log2(TPMA / TPMB)',
       x = 'Species with oCGI and Peak') +
  coord_cartesian(ylim = c(-1.9, 1.9)) +
  geom_text(data = RNAtable %>%
              filter(speciesPair %in% c('Rhesus\nRat',
                                        'Mouse\nPig',
                                        'Rat\nPig',
                                        'Pig\nHorse',
                                        'Cat\nHorse'),
                     tissue == 'Brain'),
            aes(x = x , y = y, label = stars), inherit.aes = F, size = 6)

ggsave('Fig6_RNA.pdf', Fig6, height = 1400, width = 1600, units = 'px')


################################
# Distributions for supplement #
################################

# DISTRIBUTION OF TPM VALUES - EXAMPLE IS RAT VS PIG (Fig S35B)

table$TPM_Avg <- sqrt( table$TPM.A * table$TPM.B)

table$CGI_category <- table$CGI_summary
table[is.na(table$CGI_summary),]$CGI_category <- 'Background'
table[table$CGI_category=='A',]$CGI_category <- 'A-only'
table[table$CGI_category=='B',]$CGI_category <- 'B-only'
x <- table %>% dplyr::filter(group == 'rn7_susScr11_brain_H3K27ac',
                             CGI_category != 'mix')

library(ggridges)

x$CGI_category <- factor(x$CGI_category,
                         levels = rev(c('A-only', 'B-only', 'Background')))

supplement_distribution <- x %>%
  ggplot(aes(x = log10(TPM_Avg), y = CGI_category, fill = CGI_summary)) +
  stat_density_ridges(quantile_lines = T, quantiles = 2) +
  scale_fill_manual(values = c('gray30', 'lightskyblue', 'royalblue4')) +
  theme_bw() +
  theme(legend.position = 'none')

ggsave('FigS35_TPM_distribution.pdf', supplement_distribution, height = 1000, width = 1500, units = 'px')
  

# DISTRIBUTION OF RESAMPLING MEDIANS (Fig S35D)

resamplingMedians_A <- read.table('~/Desktop/CGI/RNA/rn7_susScr11_brain_H3K27ac_A_resamplingMedians.txt', header = F) %>% tibble()
resamplingMedians_B <- read.table('~/Desktop/CGI/RNA/rn7_susScr11_brain_H3K27ac_B_resamplingMedians.txt', header = F) %>% tibble()
observedMedians <- read.table('~/Desktop/CGI/RNA/rn7_susScr11_brain_H3K27ac_observedMedians.txt', header = F) %>% tibble()

observedMedian_A <- observedMedians[,1]
observedMedian_B <- observedMedians[,2]

pdf(file = 'FigS35_exampleHistogram_A.pdf', height = 3, width = 4)
hist(resamplingMedians_A$V1,
     xlim = c(0.8, 2.0),
     las = 1,
     xlab = '',
     main = '')
abline(v = median(resamplingMedians_A$V1), col = 'black', lwd = 2)
abline(v = observedMedian_A, col = 'royalblue4', lwd = 2)
dev.off()

pdf(file = 'FigS35_exampleHistogram_B.pdf', height = 3, width = 4)
hist(resamplingMedians_B$V1,
     xlim = c(0.8, 2.0),
     las = 1,
     xlab = '',
     main = '')
abline(v = median(resamplingMedians_B$V1), col = 'black', lwd = 2)
abline(v = observedMedian_B, col = 'lightskyblue', lwd = 2)
dev.off()

# get exact normalized log2(TPM ratio) for the text
View(RNAtable %>% filter(speciesPair == 'Rat\nPig',
                         mark == 'H3K27ac',
                         tissue == 'Brain'))
log2(1.8143026/1.386228) # A-only: 0.3882505 --> 0.39
log2(0.8975093/1.396635) # B-only: -0.6379562 --> -0.64
