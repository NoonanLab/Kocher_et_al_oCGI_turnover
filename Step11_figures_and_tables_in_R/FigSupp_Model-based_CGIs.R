# 3/5/2024
# Purpose: make plots exploring relationship between my UCSC_AL CGIs, standard UCSC CGIs, and Model-based CGIs

# MAKES REVISED FIG S2

# Import libraries
library(tidyverse)
library(ggplot2)
library(cowplot)
library(patchwork)
theme_set(theme_cowplot())

# Load directory and data
setwd('/Users/acadiak/Desktop/Yale/!Writing/Kocher_CpG/GB_Revisions/Revised_code/Step12_revisions/Model-based_CGIs/')
filePathToData <- '/Users/acadiak/Desktop/Yale/!Writing/Kocher_CpG/GB_Revisions/Revised_code/Step12_revisions/Model-based_CGIs/CGI_comparisons_20240319'

# store file names in a list
#files <- dir(path = filePathToData, pattern = "*.txt")
files <- list.files(path = filePathToData,
                    pattern = "*.txt",
                    full.names = T)

files2 <- data.frame(fileName = str_replace(files, paste0(filePathToData, '/'), ''), 
                     id = as.character(1:length(files)))
# import data
data <- files %>% 
  lapply(read_tsv) %>%
  bind_rows(.id = "id") %>%
  left_join(files2) %>%
  mutate(Species = str_replace(fileName, '_CGI_comparisons.txt', '')) %>%
  select(-id, -fileName)

data <- data %>%
  filter(Species %in% c('hg19', 'rheMac10', 'mm39', 'rn7', 'canFam6', 'equCab3'))

data$Species <- factor(data$Species,
                       levels = rev(c('hg19', 'rheMac10', 'mm39', 'rn7', 'canFam6', 'equCab3')),
                       labels = rev(c('Human', 'Rhesus', 'Mouse', 'Rat', 'Dog', 'Horse')))

# add column describing situation for all pairs
data$Model_based_vs_UCSC_AL <- paste0(data$Model_based, data$UCSC_AL)
data$Model_based_vs_UCSC_Standard <- paste0(data$Model_based, data$UCSC_Standard)
data$UCSC_AL_vs_UCSC_Standard <- paste0(data$UCSC_AL, data$UCSC_Standard)

# refactor
data$Model_based_vs_UCSC_AL <- factor(data$Model_based_vs_UCSC_AL,
                                      levels = c('10', '11', '01'),
                                      labels = c('Only HMM',
                                                 'Both',
                                                 'Only UCSC permissive'))

data$Model_based_vs_UCSC_Standard <- factor(data$Model_based_vs_UCSC_Standard,
                                            levels = c('10', '11', '01'),
                                            labels = c('Only HMM',
                                                       'Both',
                                                       'Only UCSC standard'))

data$UCSC_AL_vs_UCSC_Standard <- factor(data$UCSC_AL_vs_UCSC_Standard,
                                            levels = c('10', '11', '01'),
                                            labels = c('Only UCSC permissive',
                                                       'Both',
                                                       'Only UCSC standard'))

# add column summarizing percent of length that is repeats
data$percentRepeats <- data$Rmsk_bp / data$Length * 100
data$percentNs <- data$N_number / data$Length * 100
  

################
#### plots #####
################

#####
# A #
#####

# UCSC_AL vs UCSC_Standard
A1 <- data %>%
  filter(UCSC_AL == 1 | UCSC_Standard == 1,
         oCGI == 1) %>%
  group_by(Species, UCSC_AL_vs_UCSC_Standard) %>%
  count() %>%
  ggplot(aes(x = Species, y = n/1000, fill = as_factor(UCSC_AL_vs_UCSC_Standard))) +
  geom_col() + 
  ylim(c(0, 100)) + 
  scale_fill_manual(values = c('dodgerblue4', 'gray60', 'firebrick3')) +
  theme_bw() +
  theme(legend.position="bottom",
        panel.grid.minor.y = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        strip.background = element_blank(),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 12)) +
  coord_flip() +
  labs(y = 'Number of oCGIs',
       x = '',
       fill = '')

# Model-based vs UCSC_Standard
A2 <- data %>%
  filter(Model_based == 1 | UCSC_Standard == 1,
         oCGI == 1) %>%
  group_by(Species, Model_based_vs_UCSC_Standard) %>%
  count() %>%
  ggplot(aes(x = Species, y = n/1000, fill = as_factor(Model_based_vs_UCSC_Standard))) +
  geom_col() + 
  ylim(c(0, 100)) +
  scale_fill_manual(values = c('mediumseagreen', 'gray60', 'firebrick3')) +
  theme_bw() +
  theme(legend.position="bottom",
        panel.grid.minor.y = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        strip.background = element_blank(),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 12)) +
  coord_flip() +
  labs(y = 'Number of oCGIs',
       x = '',
       fill = '')

#####
# B #
#####

# Model-based vs UCSC_AL
B <- data %>%
  filter(Model_based == 1 | UCSC_AL == 1,
         oCGI == 1) %>%
  group_by(Species, Model_based_vs_UCSC_AL) %>%
  count() %>%
  ggplot(aes(x = Species, y = n/1000, fill = as_factor(Model_based_vs_UCSC_AL))) +
  geom_col() + 
  ylim(c(0, 125)) +
  scale_fill_manual(values = c('mediumseagreen', 'gray60', 'dodgerblue4')) +
  theme_bw() +
  theme(legend.position="bottom",
        panel.grid.minor.y = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        strip.background = element_blank(),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 12)) +
  coord_flip() +
  labs(y = 'Number of oCGIs',
       x = '',
       fill = '')

#####
# C #
#####

# plot repeat content
C <- data %>%
  filter(Model_based == 1 | UCSC_AL == 1,
         oCGI == 1,
         Species %in% c('Human', 'Mouse', 'Horse')) %>%
  ggplot(aes(x = percentNs, fill = Model_based_vs_UCSC_AL)) +
  geom_histogram(binwidth = 2) +
  facet_grid(Model_based_vs_UCSC_AL ~ factor(Species, levels = c('Human', 'Mouse', 'Horse')), scales = 'free') +
  theme_bw() +
  scale_fill_manual(values = c('seagreen', 'gray60', 'dodgerblue4')) +
  theme(legend.position="bottom",
        panel.grid.minor.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(size = 12),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 12)) +
  labs(y = 'Number of oCGIs',
       x = 'Percent N bases',
       fill = '')

#####
# D #
#####

# Now filter sites with a lot of repeat content and compare HMM vs UCSC permissive
D <- data %>%
  filter(Model_based == 1 | UCSC_AL == 1,
         oCGI == 1,
         percentNs < 50) %>%
  group_by(Species, Model_based_vs_UCSC_AL) %>%
  count() %>%
  ggplot(aes(x = Species, y = n/1000, fill = as_factor(Model_based_vs_UCSC_AL))) +
  geom_col() + 
  ylim(c(0, 125)) +
  scale_fill_manual(values = c('mediumseagreen', 'gray60', 'dodgerblue4')) +
  theme_bw() +
  theme(legend.position="bottom",
        panel.grid.minor.y = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        strip.background = element_blank(),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 12)) +
  coord_flip() +
  labs(y = 'Number of oCGIs',
       x = '',
       fill = '')

ABD <- A1 + A2 + B + D + plot_layout(ncol = 4)

ggsave('FigS2_oCGI_overlap_ABD.pdf', ABD, height = 1400, width = 3200, units = 'px')
ggsave('FigS2_oCGI_overlap_C.pdf', C, height = 1800, width = 1800, units = 'px')


