# 2/27/2024
# Purpose: make stacked bar plots showing distribution of oCGIs across cCRE types

# MAKES FIG S6

# Import libraries
library(cowplot)
theme_set(theme_cowplot())
library(ggplot2)
library(tidyverse)

# Load directory and data
setwd('/Users/acadiak/Desktop/Yale/!Writing/Kocher_CpG/GB_Revisions/Revised_code/Step12_revisions/cCREs')
filePathToData <- '/Users/acadiak/Desktop/Yale/!Writing/Kocher_CpG/GB_Revisions/Revised_code/Step12_revisions/cCREs/oCGI_cCRE_overlap.txt'
data <- read.table(filePathToData, header = T) %>%
  tibble()

# Add percent column
data <- data %>%
  group_by(Species, Feature) %>%
  mutate(percent =  n/sum(n)) %>% 
  ungroup

# Store order of factors
Species <- c('Mouse', 'Human')
cCREs <- c('PLS', 'PLS,CTCF-bound', 'pELS', 'pELS,CTCF-bound', 'dELS', 'dELS,CTCF-bound', 'DNase-H3K4me3', 'DNase-H3K4me3,CTCF-bound', 'CTCF-only,CTCF-bound')
Features <- c('All', 'Intronic-intergenic', 'oCGIs')

# Reorder cCRE levels
data$cCRE <- factor(data$cCRE,
       levels = cCREs,
       labels = c('PLS', 'PLS + CTCF', 'pELS', 'pELS + CTCF', 'dELS', 'dELS + CTCF', 'DNase-H3K4me3', 'DNase-H3K4me3 + CTCF', 'CTCF only'))

# Relabel "features"
data$Feature <- factor(data$Feature,
                    levels = Features,
                    labels = c('All cCREs', 'Intronic-intergenic cCREs', 'cCREs overlapping oCGIs'))

# Plot
cCRE_overlap <- data %>%
  ggplot(aes(x = factor(Feature), y = percent, fill = cCRE)) +
  geom_col() +
  facet_grid(. ~ Species) +
  theme_bw() +
  theme(panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        strip.background = element_blank(),
        axis.text = element_text(size = 12),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
        axis.title = element_text(size = 12),
        strip.text = element_text(size = 12)) +
  scale_y_continuous(labels = scales::percent,
                     breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1.0),
                     expand = c(0.01, 0.01),
                     limits = c(0, 1.0)) +
  scale_fill_manual(values = c('firebrick2', 'darkred', 'darkslategray2', 'darkslategray3', 'darkseagreen1', 'darkseagreen', 'orange', 'darkorange2', 'dodgerblue3')) +
  labs(y = 'Percent',
       x = '')

ggsave('FigS6_cCRE_overlap.pdf', cCRE_overlap, height = 1800, width = 2800, units = 'px')

