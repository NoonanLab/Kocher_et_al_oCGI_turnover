# 10/21/22
# Purpose: make singleSpecies CGI-centric parts of Fig 1 (1A, 1B)
# using summary tables made in Step3
# updated 1/2/23 to update figure formatting and include permutation test for 1B

# load required packages
require(ggplot2)
require(tidyverse)
require(cowplot)
require(purrr)
theme_set(theme_cowplot())
library(patchwork)

# store directory to save figures and file path to data tables (generated in Step3)
setwd('/Users/acadiak/Desktop/CGI/Figures/Mar_2023')
filePathToData <- '/Users/acadiak/Desktop/CGI/singleSpecies/summaryTables_singleSpeciesCGIcentric/'

# store file names in a list
files <- dir(path = filePathToData, pattern = "*.txt")

# read in as nested dataframe
data <- tibble(filename = files) %>%
  mutate(file_contents = map(filename,
                             function(x) read_tsv(file.path(filePathToData, x)))
  )

# unnest to generate full table, but with file name as a column
unnested_CGI <- unnest(data, cols = c(file_contents))

# modify variable names using various tidy commands
preFilter_CGI <- unnested_CGI %>% 
  # modify file names labels and order factors
  mutate(filename = factor(filename, 
                            levels = files, 
                            labels = c('Marmoset', 'Dog', 'Horse', 'Cat', 'Human', 'Mouse', 'Rhesus', 'Rat', 'Pig'))) %>%
  # modify mark labels and order factors
  mutate_at(.vars = vars(CGI_liftsToHg38), factor) # remove this: starts_with('peak.')

# filter to only sites that lift to human and have no feature in human (both are incorporated in variable CGI_liftsToHg38)
table_CGI <- preFilter_CGI %>%
  filter(CGI_liftsToHg38=='1')

#############################################
# generate Fig 1A - counts of oCGIs in each species

# reorder species
table_CGI$species <- fct_relevel(table_CGI$filename, rev(c('Human', 'Rhesus', 'Marmoset', 
                                                           'Mouse','Rat', 'Pig', 
                                                           'Dog', 'Cat', 'Horse')))

# plot number of oCGIs in each species
table_CGI %>%
  group_by(species) %>%
  count() %>%
  ggplot(aes(species, n/1000)) +
  geom_bar(stat = 'identity', width = 0.6) + 
  coord_flip() + 
  theme_minimal() +
  labs(y = 'Number of oCGIs (in Thousands)',
       x = '') +
  theme(axis.line = element_line(color = 'black', size = 1),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_line(size = 1),
        text = element_text(size = 20),
        axis.text = element_text(size = 36)) +
  scale_y_continuous(limits = c(0, 85), expand = c(0, 0))

ggsave('Fig1A_CGIcounts_barplot.pdf', height = 2800, width = 2150, units = 'px')

# view counts in each species
View(table_CGI %>%
  group_by(species) %>%
  count())
# mouse has fewest (11067) and cat has most (77199)
# Mouse 11067
# Marmoset 18238
# Rhesus 21848
# Rat 23717
# Human 24487
# Dog 37905
# Pig 49737
# Horse 53260
# Cat 77199

###############################################################
# generate Fig 1B - counts with activity in each species

shuffling_CGIcentric <- read_tsv('/Users/acadiak/Desktop/CGI/singleSpecies/reshufflingSummary_CGI-centric.txt',
                                 col_names = F) %>% tibble()
colnames(shuffling_CGIcentric) <- c('species', 'tissue', 'mark', 'observed', 'expected', 'p')

# add column with CGI vs peak called "centric"
shuffling_CGIcentric$centric <- 'CGI'

# add line with extra species to make plot alignment work
shuffling_CGIcentric <- shuffling_CGIcentric %>%
  add_row(species = 'hg19', tissue = 'brain', mark = 'H3K4me3',
          observed = 0.0, expected = 0.0, p = 1.0, centric = 'CGI') %>%
  add_row(species = 'hg19', tissue = 'brain', mark = 'H3K4me1',
          observed = 0.0, expected = 0.0, p = 1.0, centric = 'CGI')
for (species in c('calJac4', 'rn7', 'susScr11', 'canFam6', 'felCat9', 'equCab3')) {
  shuffling_CGIcentric <- shuffling_CGIcentric %>%
    add_row(species = species, tissue = 'devBrain', mark = 'H3K4me2',
            observed = 0.0, expected = 0.0, p = 1.0, centric = 'CGI')
}

# rename table
shuffling <- shuffling_CGIcentric

# adjust p values
shuffling$padj <- p.adjust(shuffling$p, method = 'BH')

# rename species, tissue, and mark
shuffling <- shuffling %>%
  mutate(species = factor(species,
                          levels = c('hg19', 'rheMac10', 'calJac4', 'mm39', 'rn7', 'susScr11', 'canFam6', 'felCat9', 'equCab3'),
                          labels = c('Human', 'Rhesus', 'Marmoset', 'Mouse', 'Rat', 'Pig', 'Dog', 'Cat', 'Horse')),
         tissue = factor(tissue,
                       levels = c('brain', 'liver', 'muscle', 'testis', 'devBrain', 'devLimb', 'all'),
                       labels = c('B', 'L', 'M', 'T', 'DB', 'DL', 'Any')),
         mark = factor(mark,
                       levels = c('H3K4me3', 'H3K27ac', 'H3K4me1', 'H3K4me2', 'all'),
                       labels = c('H3K4me3', 'H3K27ac', 'H3K4me1', 'H3K4me2', 'Any')))

# factor variables to get plots in the right order
shuffling$species <- factor(shuffling$species, levels = c('Human', 'Marmoset', 'Mouse', 'Rat', 'Pig', 'Dog', 'Cat', 'Horse', 'Rhesus'))
shuffling$tissue <- factor(shuffling$tissue, levels = c('B', 'L', 'M', 'T', 'DB', 'DL', 'Any'))
shuffling$mark <- factor(shuffling$mark, levels = c('H3K4me3', 'H3K27ac', 'H3K4me1', 'H3K4me2', 'Any'))

# add values for geom_segment
shuffling$x = 0
shuffling[shuffling$tissue == 'B',]$x <- 0.5
shuffling[shuffling$tissue == 'L',]$x <- 1.5
shuffling[shuffling$tissue == 'M',]$x <- 2.5
shuffling[shuffling$tissue == 'T',]$x <- 3.5
shuffling[shuffling$tissue == 'DB',]$x <- 4.5
shuffling[shuffling$tissue == 'DL',]$x <- 5.5
shuffling[shuffling$tissue == 'Any',]$x <- 0.5
shuffling[shuffling$tissue == 'DB' & shuffling$mark == 'H3K4me2',]$x <- 0.5
shuffling$xend <- shuffling$x + 1
shuffling$y <- shuffling$expected
shuffling$yend <- shuffling$y
shuffling$stars <- NA
shuffling[shuffling$padj<0.01,]$stars <- '*'
shuffling$stars.x <- (shuffling$xend + shuffling$x) / 2

# define function to make barplots
cp.marks <- c('darkorange1', 'palegreen4', 'mediumpurple4', 'darkgoldenrod1', 'gray30')

function_barplots <- function(table, inMark, color, ylimit) {
  x <- table %>%
    filter(mark == inMark) %>%
    ggplot() +
    facet_grid(species ~ ., scales = 'free', space = 'free') +
    geom_col(aes(x = tissue, y = observed, fill = color), width = 0.8) +
    scale_fill_manual(values = color, guide = 'none') +
    theme_bw() +
    theme(panel.grid.minor.x = element_blank(),
          panel.grid.major.x = element_blank(),
          panel.grid.minor.y = element_blank(),
          strip.background = element_blank(),
          axis.text = element_text(size = 12),
          strip.text = element_text(size = 16)) +
    scale_y_continuous(labels = scales::percent,
                       limits = c(0, ylimit*1.1),
                       breaks = c(0, ylimit / 2, ylimit),
                       expand = c(0, 0)) +
    labs(y = '',
         x = '') +
    geom_segment(data = table %>% filter(mark == inMark,
                                         stars == '*'),
                 aes(x = x, xend = xend,
                     y = y, yend = yend), 
                   inherit.aes = F,
                   lwd = 0.5,
                 color = 'gray80') +
    geom_text(data = table %>% filter(mark == inMark,
                                      stars == '*'),
              aes(x = stars.x , y = observed + 0.02, label = stars), # y = ylimit * 0.75
              inherit.aes = F,
              size = 6)
  return(x)
}

CGI_A <- function_barplots(shuffling %>% filter(centric == 'CGI', species != 'Rhesus'), 'H3K4me3', cp.marks[1], 0.2)
CGI_B <- function_barplots(shuffling %>% filter(centric == 'CGI', species != 'Rhesus'), 'H3K27ac', cp.marks[2], 0.50)
CGI_C <- function_barplots(shuffling %>% filter(centric == 'CGI', species != 'Rhesus'), 'H3K4me1', cp.marks[3], 0.50)
CGI_D <- function_barplots(shuffling %>% filter(centric == 'CGI', species != 'Rhesus'), 'H3K4me2', cp.marks[4], 0.50)
CGI_E <- function_barplots(shuffling %>% filter(centric == 'CGI', species != 'Rhesus'), 'Any', cp.marks[5], 0.80)

# generate Fig S5 (originally S3)
CGI_centric <- CGI_A + CGI_B + CGI_C + CGI_D + CGI_E + plot_layout(ncol = 5, widths = c(4, 6, 4, 1, 1))
ggsave('FigS3_allSpecies_CGIcentric.pdf', CGI_centric, height = 2700, width = 2750, units = 'px')

# get values for text
View(shuffling %>% filter(tissue == 'Any', mark == 'Any'))

##########
# Fig 1B #
##########

# make rhesus only with different proportions for main figure
Rhesus_CGI_A <- function_barplots(shuffling %>% filter(centric == 'CGI', species == 'Rhesus'), 'H3K4me3', cp.marks[1], 0.2)
Rhesus_CGI_B <- function_barplots(shuffling %>% filter(centric == 'CGI', species == 'Rhesus'), 'H3K27ac', cp.marks[2], 0.50)
Rhesus_CGI_C <- function_barplots(shuffling %>% filter(centric == 'CGI', species == 'Rhesus'), 'H3K4me1', cp.marks[3], 0.50)
Rhesus_CGI_D <- function_barplots(shuffling %>% filter(centric == 'CGI', species == 'Rhesus'), 'H3K4me2', cp.marks[4], 0.50)
Rhesus_CGI_E <- function_barplots(shuffling %>% filter(centric == 'CGI', species == 'Rhesus'), 'Any', cp.marks[5], 0.80)

# JUST DO CGI-CENTRIC IN MAIN FIGURE
Fig1B <- Rhesus_CGI_A + Rhesus_CGI_B + Rhesus_CGI_C + Rhesus_CGI_D + Rhesus_CGI_E + plot_layout(ncol = 5, widths = c(4, 6, 4, 1, 1))
ggsave('Fig1B_rhesus.pdf', Fig1B, height = 580, width = 2700, units = 'px')



