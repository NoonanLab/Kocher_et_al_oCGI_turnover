# 11/3/22
# Purpose: analyze ChIP-seq data from hs754 humanized mouse 
# using DESeq2 to identify differential peaks between WT and HUM

require(DESeq2)
require(tidyverse)
require(purrr)
require(cowplot)
theme_set(theme_cowplot())
require(scattermore)

setwd('/Users/acadiak/Desktop/CGI/hs754/ChIP')

# define function for running DESeq2, outputting results which will be assigned as a separate object for each analysis
runDESeq2 <- function(multiplexGroup, timePoint, mark) {
  fileStem <- paste(multiplexGroup, timePoint, mark, sep = '_')
  
  # store file names
  fileNames <- list.files(path = 'counts',
                          pattern = fileStem,
                          full.names = T)
  sampleNames <- sub('counts/', '', fileNames)
  sampleNames <- sub('.quant', '', sampleNames)
  
  # get genotypes
  sampleInfo <- data.frame(sampleName=sampleNames,
                           fileName=fileNames,
                           genotype=rep('WT', length(fileNames)))
  sampleInfo[str_detect(fileNames, "HUM"),]$genotype <- "HUM"
    
  # read in data
  dds <- DESeqDataSetFromHTSeqCount(sampleTable=sampleInfo,
                                    directory=getwd(),
                                    design=~genotype)
  dds2 <- DESeq(dds)
  res <- results(dds2, contrast = c('genotype', 'HUM', 'WT'))
  return(res)
}

# run DESeq2 for all combinations
res_e17.5_H3K27ac <- runDESeq2('D', 'e17.5', 'H3K27ac')
res_e17.5_H3K4me3 <- runDESeq2('A', 'e17.5', 'H3K4me3')
res_e17.5_CTCF <- runDESeq2('A', 'e17.5', 'CTCF')
res_e11.5_H3K27ac <- runDESeq2('E', 'e11.5', 'H3K27ac')
res_e11.5_H3K4me3 <- runDESeq2('B', 'e11.5', 'H3K4me3')

# convert to data frames and add columns with name info
res_e17.5_H3K27ac$label <- 'e17.5_H3K27ac'
res_e17.5_H3K4me3$label <- 'e17.5_H3K4me3'
res_e17.5_CTCF$label <- 'e17.5_CTCF'
res_e11.5_H3K27ac$label <- 'e11.5_H3K27ac'
res_e11.5_H3K4me3$label <- 'e11.5_H3K4me3'

# paste together results from each separate DESeq run
fullTable <- bind_rows(data.frame(res_e17.5_H3K27ac),
                       data.frame(res_e17.5_H3K4me3),
                       data.frame(res_e17.5_CTCF),
                       data.frame(res_e11.5_H3K27ac),
                       data.frame(res_e11.5_H3K4me3),
                       ) %>%
  separate(label, c('timePoint', 'mark'), sep = '_') %>%
  rownames_to_column(var = 'coord') %>%
  separate(coord, c('chr', 'temp'), sep = ':', remove = F) %>%
  separate(temp, c('start', 'end'), sep = '-') %>%
  drop_na()

# specify alpha for plotting - very transparent, but opaque if significant genome-wide
fullTable$alpha <- 0.05
fullTable[fullTable$padj<0.05,]$alpha <- 1

# specify category for changing color:
# background, genome-wide significant (any, on chr19, and in hs754), and hs754 significant but not genome-wide
fullTable$category <- 'background'
fullTable[fullTable$padj<0.05,]$category <- 'genome-wide significant'
fullTable[fullTable$padj<0.05&grepl('chr19', fullTable$chr),]$category <- 'genome-wide significant and on chr19'
fullTable$hs754 <- 0
fullTable[fullTable$chr=='chr13'&fullTable$start>72436073&fullTable$start<72441214,]$hs754 <- 1
fullTable[fullTable$hs754==1&fullTable$pvalue<0.05&fullTable$padj>0.05,]$category <- 'hs754, significant but not genome-wide'
fullTable[fullTable$hs754==1&fullTable$pvalue<0.05&fullTable$padj>0.05,]$alpha <- 1
fullTable[fullTable$hs754==1&fullTable$padj<0.05,]$category <- 'genome-wide significant and hs754'

# specify factor order
fullTable$timePoint <- factor(fullTable$timePoint, levels = c('e11.5', 'e17.5'))
fullTable$mark <- factor(fullTable$mark, levels = c('H3K4me3', 'H3K27ac', 'H3K27me3', 'CTCF'))

# view peaks in hs754
View(fullTable[fullTable$chr=='chr13'&fullTable$start>72436073&fullTable$start<72441214,])

# make plots for Fig S32
ChIPplots <- fullTable %>%
  filter(mark %in% c('H3K4me3', 'H3K27ac')) %>%
  ggplot(aes(x = log2(baseMean), 
             y = log2FoldChange, 
             color = category)) + 
  facet_grid(mark ~ timePoint) +
  coord_cartesian(ylim = c(-3.5, 3.5)) +
  geom_hline(yintercept = 0) +
  geom_scattermore(aes(alpha = alpha),
                   pointsize = 10.2,
                   pixels = c(1200, 1200)) + 
  scale_color_manual(values = c('gray30', 'orange1','firebrick3', 'gold1', 'hotpink')) +
  theme_bw() +
  theme(panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        strip.background = element_blank(),
        axis.text = element_text(size = 20),
        strip.text = element_text(size = 20),
        axis.title = element_text(size=20)) +
  labs(x = 'Log2 (baseMean)',
       y = 'Log2 (HUM / WT)')

setwd('/Users/acadiak/Desktop/CGI/Figures/Jan_2023')
ggsave('FigS32_hs754_ChIP.pdf', ChIPplots, height = 2400, width = 3300, units = 'px')

# get table for significantly differential peaks
View(fullTable %>% filter(padj < 0.05))

###################
# export Table S9 #
###################

# change category column to just be NA vs 'hs754' vs 'chr19 CNV'
fullTable$category <- NA
fullTable[fullTable$chr=='chr19'&((fullTable$start>37215255&fullTable$start<37327255)|(fullTable$start>36889256&fullTable$start<37369503)),]$category <- 'chr19 CNV'
fullTable[fullTable$chr=='chr13'&fullTable$start>72436073&fullTable$start<72441214,]$category <- 'hs754'

# column with whether padj < 0.05
fullTable$significant <- 'no'
fullTable[fullTable$padj<0.05,]$significant <- 'yes'

# move columns around and rename some
fullTableForExport <- fullTable %>%
  filter(mark %in% c('H3K27ac', 'H3K4me3', 'CTCF')) %>%
  select(c('timePoint', 'mark', 'coord', 'category', 'baseMean', 'log2FoldChange', 'pvalue', 'padj', 'significant')) %>%
  mutate(timePoint = factor(timePoint,
                            levels = c('e11.5', 'e17.5'),
                            labels = c('E11.5', 'E17.5')),
         mark = factor(mark,
                       levels = c('H3K27ac', 'H3K4me3', 'CTCF')),
         significant = factor(significant,
                              levels = c('yes', 'no')))

# rename columns
colnames(fullTableForExport) <- c('TimePoint', 'HistoneModification_orTranscriptionFactor',
                                  'PeakCoord', 'Description', 'baseMean Counts',
                                  'Log2(Fold Change)', 'pvalue', 'qvalue', 'Significant')

# write to table
setwd('/Users/acadiak/Desktop/CGI/Figures/Mar_2023')
write_delim(fullTableForExport %>% arrange(Significant, TimePoint, HistoneModification_orTranscriptionFactor,
                                           qvalue, pvalue, PeakCoord),
            file = 'TableS9_ChIP_results.txt',
            delim = '\t',
            col_names = T)


# export subset of table as bed files for the chr19 CNV region (for Fig S43)
# viewing region: chr19:36,748,165-37,454,636

View(fullTable[fullTable$chr=='chr13'&fullTable$start>72436073&fullTable$start<72441214,])
# filter to region I want to include in the figure
bedTable <- fullTable %>%
  filter(chr == 'chr19',
         start > 36748165,
         start < 37454636)

# make new category column showing significance and genome-wide significance (currently only have label for genome-wide)
# REDO on 5/4/23 with p < 0.05 (instead of 0.01 as done previously) - use to recolor peaks in illustrator file (only turned a few gray peaks yellow)
bedTable$category2 <- bedTable$category
bedTable[bedTable$pvalue < 0.05,]$category2 <- 'significant'
bedTable[bedTable$padj < 0.05,]$category2 <- 'genome-wide significant'

# add new columns for all bed columns
# chr start end log2FC score(0) strand(.) thickStart(start) thickEnd(end) RGB
# RGB for background = gray rgb(140, 140, 140)
# RGB for significant = gold1 rgb(255,215,0)
# RGB for genome-wide significant = darker orange rgb(255, 101, 23)
bedTable$score <- 0
bedTable$strand <- '.'
bedTable$thickStart <- bedTable$start
bedTable$thickEnd <- bedTable$end
bedTable$rgb <- '140,140,140'
bedTable[bedTable$category2=='significant',]$rgb <- '255,215,0'
bedTable[bedTable$category2=='genome-wide significant',]$rgb <- '255,101,23'

# write bed files
setwd('~/Desktop/CGI/hs754/CNV')
for (tp in c('e11.5', 'e17.5')) {
  for (mk in c('H3K4me3', 'H3K27ac')) {
    bedTable %>%
      filter(timePoint == tp,
             mark == mk) %>%
      select(chr, start, end, log2FoldChange, score, strand, thickStart, thickEnd, rgb) %>%
      write_delim(file = paste0(tp, '_', mk, '_CNVregion.bed'),
                  delim = '\t',
                  col_names = F)
  }
}

# use in browser to generate Figure S5.3
# before each bed file paste this line:
# track type=bed name=e11.5_H3K4me3 itemRgb=On
# track type=bed name=e11.5_H3K27ac itemRgb=On
# track type=bed name=e17.5_H3K4me3 itemRgb=On
# track type=bed name=e17.5_H3K27ac itemRgb=On

