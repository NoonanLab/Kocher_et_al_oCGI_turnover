# 10/27/22
# Run DESeq2 on RNA-seq data from developing diencephalon
# 4 replicates per genotype (WT and HUM) at both e11.5 and e17.5

# THIS WAS FORMERLY ASSOCIATED WITH FIGURE 5 - NOW REVISED FIGURE 4

#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("DESeq2")

require(DESeq2)
require(tidyverse)
require(purrr)
require(patchwork)
require(scattermore)

setwd('/Users/acadiak/Desktop/CGI/hs754/RNA')

# RUN FOR e11.5
fileNames <- list.files(path = 'results_hs754_RNA-seq',
                        pattern = "e11.5",
                        full.names = T)

# load data into wideTable where each col is a sample and each row is a gene
wideTable <- fileNames %>% 
  map_dfr(read_tsv, 
          col_names = c('geneID', 'counts.unstranded', 'counts.forwardStrand', 'counts.reverseStrand'),
          col_select = c('geneID', 'counts.reverseStrand'),
          id = 'file') %>%
  pivot_wider(names_from = 'file',
              values_from = 'counts.reverseStrand') %>%
  # turn geneID column into rownames
  column_to_rownames(var = 'geneID') %>%
  # remove first 4 rows
  filter(!row_number() %in% c(1,2,3,4))

colnames(wideTable)

# store sampleInfo in same order as columns in wideTable
sampleInfo <- data.frame(genotype=c(rep('HUM', 4),
                                    rep('WT', 4)))

# run DESeq2
dds11 <- DESeqDataSetFromMatrix(countData = wideTable,
                              colData = sampleInfo, 
                              design = ~ genotype)

dds11 <- DESeq(dds11)
res11 <- results(dds11, contrast = c('genotype', 'HUM', 'WT'), alpha = 0.05)

# view results
summary(res11)
head(res11[order(res11$padj),])

# convert to data frame and add columns with alpha and color
# drop rows with log2FoldChange or padj set to NA
# see here for info on why that might be http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#pvaluesNA
res11 <- data.frame(res11) %>% drop_na(log2FoldChange, padj)
res11$a <- 0.05
res11[res11$padj<0.05,]$a <- 1

# modify colors for IRX1 and IRX2
res11['ENSMUSG00000060969',]$a <- 0.98
res11['ENSMUSG00000001504',]$a <- 0.98

# in ggplot
e11.5 <- res11 %>%
  ggplot(aes(x = log2(baseMean), y = log2FoldChange, alpha = a, color = factor(a))) +
  geom_scattermore(aes(alpha = a),
                   pointsize = 7,
                   pixels = c(1200, 1200)) +
  theme_bw() +
  scale_color_manual(values = c('gray30', 'deepskyblue3', 'deeppink')) +
  theme(axis.line = element_line(color = 'black', size = 1),
        panel.grid.minor.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        text = element_text(size = 16),
        axis.text = element_text(size = 16),
        legend.position = 'none') + 
  labs(x = 'Log2(Average Counts)',
       y = 'Log2(HUM / WT)') +
  coord_cartesian(ylim = c(-5, 5))
# top significant gene is ENSMUSG00000022092 = Ppp3cc on chr 14
# bottom significant gene is ENSMUSG00000021703 = Serinc5 on chr 13, approx 20 Mb from hs754
# coord = chr13:92747599-92848455
# distance from hs754 = 92747599 - 72441214 = 20.3 Mb

# left blue dot is Irx1, right blue dot is Irx2


# RUN FOR e17.5
fileNames <- list.files(path = 'results_hs754_RNA-seq',
                        pattern = "e17.5",
                        full.names = T)

# load data into wideTable where each col is a sample and each row is a gene
wideTable <- fileNames %>% 
  map_dfr(read_tsv, 
          col_names = c('geneID', 'counts.unstranded', 'counts.forwardStrand', 'counts.reverseStrand'),
          col_select = c('geneID', 'counts.reverseStrand'),
          id = 'file') %>%
  pivot_wider(names_from = 'file',
              values_from = 'counts.reverseStrand') %>%
  # turn geneID column into rownames
  column_to_rownames(var = 'geneID') %>%
  # remove first 4 rows
  filter(!row_number() %in% c(1,2,3,4))

colnames(wideTable)

# store sampleInfo in same order as columns in wideTable
sampleInfo <- data.frame(genotype=c(rep('HUM', 4),
                                    rep('WT', 4)))

# run DESeq2
dds17 <- DESeqDataSetFromMatrix(countData = wideTable,
                              colData = sampleInfo, 
                              design = ~ genotype)

dds17 <- DESeq(dds17)
res17 <- results(dds17, contrast = c('genotype', 'HUM', 'WT'), alpha = 0.05)

# view results
summary(res17)
res17[order(res17$padj),]

# convert to data frame and add columns with alpha and color
# drop rows with log2FoldChange or padj set to NA
# see here for info on why that might be http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#pvaluesNA
res17 <- data.frame(res17) %>% drop_na(log2FoldChange, padj)
res17$a <- 0.05
res17[res17$padj<0.05,]$a <- 1

# modify colors for IRX1 and IRX2
res17['ENSMUSG00000060969',]$a <- 0.98
res17['ENSMUSG00000001504',]$a <- 0.98

# in ggplot
e17.5 <- res17 %>%
  ggplot(aes(x = log2(baseMean), y = log2FoldChange, alpha = a, color = factor(a))) +
  geom_scattermore(aes(alpha = a),
                   pointsize = 7,
                   pixels = c(1200, 1200)) +
  theme_bw() +
  scale_color_manual(values = c('gray30', 'deepskyblue3', 'deeppink')) +
  theme(axis.line = element_line(color = 'black', size = 1),
        panel.grid.minor.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        text = element_text(size = 16),
        axis.text = element_text(size = 16),
        legend.position = 'none') +
  labs(x = 'Log2(Average Counts)',
       y = 'Log2(HUM / WT)') +
  coord_cartesian(ylim = c(-5, 5))
# left blue dot is Irx1, right blue dot is Irx2

# print plots - Figure S40 (originally S33)
setwd('/Users/acadiak/Desktop/CGI/Figures/Jan_2023')
RNA_plots <- e11.5 + e17.5
ggsave('FigS5.4_RNA_volcanoPlots.pdf', RNA_plots, height = 1200, width = 2400, units = 'px')


# view significant DE genes for table
res11 %>% filter(padj < 0.05)
res17 %>% filter(padj < 0.05)

####################
# output Table S10 #
####################

# add time point column to both tables
res11$timePoint <- 'E11.5'
res17$timePoint <- 'E17.5'

# put gene names in a column
res11$geneName <- rownames(res11)
res17$geneName <- rownames(res17)

# combine tables and add row with significance
RNAtable <- bind_rows(res11, res17) %>% tibble()
RNAtable$significant <- NA
RNAtable[RNAtable$padj<0.05,]$significant <- 'Yes'

# manipulate for plotting
RNAtableForExport <- RNAtable %>%
  select(timePoint, geneName, baseMean, log2FoldChange, pvalue, padj, significant) %>%
  mutate(timePoint = factor(timePoint,
                            levels = c('E11.5', 'E17.5')),
         significant = factor(significant,
                              levels = c('Yes', NA)))
  
colnames(RNAtableForExport) <- c('TimePoint', 'GeneName', 'baseMean Counts',
                                 'Log2(Fold Change)', 'p-value', 'qvalue', 'Significant')

# print
setwd('/Users/acadiak/Desktop/CGI/Figures/Mar_2023')
write_delim(RNAtableForExport %>% arrange(Significant, TimePoint, qvalue, GeneName),
            file = 'TableS10_RNA_results.txt',
            delim = '\t',
            col_names = T)


