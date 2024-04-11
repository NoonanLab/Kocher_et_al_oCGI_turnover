# 3/2/2024
# Purpose: make stacked bar plots showing distribution of oCGIs across cCRE types

# MAKES REVISED FIG S21

# Import libraries
library(tidyverse)
library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())

# Load directory and data
setwd('/Users/acadiak/Desktop/Yale/!Writing/Kocher_CpG/GB_Revisions/Revised_code/Step12_revisions/Gain_vs_Loss')
filePathToData <- '/Users/acadiak/Desktop/Yale/!Writing/Kocher_CpG/GB_Revisions/Revised_code/Step12_revisions/Gain_vs_Loss/ak20240302_CGI_summary_acrossTree.txt'
data <- read.table(filePathToData, header = T) %>%
  tibble()

# Summarize species patterns as strings, with phastCons as final digit in string

data <- data %>%
  unite(col = 'pattern_all', c('hg38':'equCab3'), remove = F, sep = '')

data <- data %>%
  unite(col = 'pattern_top', c('hg38':'rn7'), remove = F, sep = '')

data <- data %>%
  unite(col = 'pattern_bottom', c('susScr11':'equCab3'), remove = F, sep = '')

### Collect numbers with certain patterns

# terminal branch-specific, with the sequence present in all
seq_all_gain <- c('211111111', '121111111', '112111111',
                  '111211111', '111121111', '111112111',
                  '111111211', '111111121', '111111112')
seq_all_loss <- c('122222222', '212222222', '221222222',
                  '222122222', '222212222', '222221222',
                  '222222122', '222222212', '222222221')

data$seq_all_gain_loss <- NA
data[data$pattern_all %in% seq_all_gain]$seq_all_gain_loss <- 'gain'

# terminal branch-specific, with the sequence present in Human / Rhesus / Marmoset / Mouse / Rat
seq_top_gain <- c('21111', '12111', '11211', '11121', '11112')
seq_top_loss <- c('12222', '21222', '22122', '22212', '22221')

# terminal branch-specific, with the sequence present in Pig / Dog / Cat / Horse
seq_bottom_gain <- c('2111', '1211', '1121', '1112')
seq_bottom_loss <- c('1222', '2122', '2212', '2221')


# Add column saying the species with gain or loss
species_all <- c('Human', 'Rhesus', 'Marmoset', 'Mouse', 'Rat', 'Pig', 'Dog', 'Cat', 'Horse')
species_top <- species_all[1:5]
species_bottom <- species_all[6:9]


################################
# even smaller species subsets #
################################

# Summarize species patterns as strings, with phastCons as final digit in string

data <- data %>%
  unite(col = 'pattern_primates', c('hg38':'calJac4'), remove = F, sep = '')

data <- data %>%
  unite(col = 'pattern_rodents', c('hg38', 'mm39', 'rn7'), remove = F, sep = '')

data <- data %>%
  unite(col = 'pattern_carnivores', c('susScr11':'felCat9'), remove = F, sep = '')

### Collect numbers with certain patterns

# terminal branch-specific within primates (gain/loss in human, gain/loss in rhesus (marmoset as outgroup))
seq_primates_gain <- c('211', '121')
seq_primates_loss <- c('122', '212')

# terminal branch-specific in rodents (gain/loss in mouse, gain/loss in rat (human as outgroup))
seq_rodents_gain <- c('121', '112')
seq_rodents_loss <- c('212', '221')

# terminal branch-specific in carnivores (gain/loss in dog, gain/loss in cat (pig as outgroup))
seq_carnivores_gain <- c('121', '112')
seq_carnivores_loss <- c('212', '221')

## make three separate tables based on whole tree, top of tree, or bottom of tree
# and add column to table saying if gain or loss

data_primates <- data %>%
  filter(pattern_primates %in% c(seq_primates_gain, seq_primates_loss)) %>%
  mutate(gain_vs_loss = case_when(
    pattern_primates %in% seq_primates_gain ~ "Gain",
    pattern_primates %in% seq_primates_loss ~ "Loss"
  ))

data_rodents <- data %>%
  filter(pattern_rodents %in% c(seq_rodents_gain, seq_rodents_loss)) %>%
  mutate(gain_vs_loss = case_when(
    pattern_rodents %in% seq_rodents_gain ~ "Gain",
    pattern_rodents %in% seq_rodents_loss ~ "Loss"
  ))

data_carnivores <- data %>%
  filter(pattern_carnivores %in% c(seq_carnivores_gain, seq_carnivores_loss)) %>%
  mutate(gain_vs_loss = case_when(
    pattern_carnivores %in% seq_carnivores_gain ~ "Gain",
    pattern_carnivores %in% seq_carnivores_loss ~ "Loss"
  ))

# Add column saying the species with gain or loss
species_primates <- species_all[1:2]
species_rodents <- species_all[4:5]
species_carnivores <- species_all[7:8]

data_primates$species <- '.'
for (i in c(1:2)) {
  data_primates[data_primates$pattern_primates==seq_primates_gain[i]|data_primates$pattern_primates==seq_primates_loss[i],]$species <- species_primates[i]
}

data_rodents$species <- '.'
for (i in c(1:2)) {
  data_rodents[data_rodents$pattern_rodents==seq_rodents_gain[i]|data_rodents$pattern_rodents==seq_rodents_loss[i],]$species <- species_rodents[i]
}

data_carnivores$species <- '.'
for (i in c(1:2)) {
  data_carnivores[data_carnivores$pattern_carnivores==seq_carnivores_gain[i]|data_carnivores$pattern_carnivores==seq_carnivores_loss[i],]$species <- species_carnivores[i]
}

# combine into single table
data_for_plots <- rbind(data_primates, data_rodents, data_carnivores)

# refactor phastCons variable
data_for_plots$phastCons <- factor(data_for_plots$phastCons,
                                  levels = c(0,1),
                                  labels = c('-', '+'))

data_for_plots$species <- factor(data_for_plots$species,
                                 levels = c('Human', 'Rhesus', 'Mouse', 'Rat', 'Dog', 'Cat'))

### Plot based on patterns

gain_loss_plot <- data_for_plots %>%
  group_by(species, gain_vs_loss, phastCons) %>%
  count() %>%
  ggplot(aes(x = as_factor(phastCons), y = n/1000, fill = gain_vs_loss)) + 
  geom_col() +
  scale_fill_manual(values = c('deepskyblue4', 'brown3')) +
  facet_grid(. ~ species, scales = 'free') +
  theme_bw() +
  theme(panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank(),
        strip.background = element_blank(),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 12),
        strip.text = element_text(size = 12)) +
  scale_y_continuous(breaks = c(0, 5, 10, 15, 20, 25),
                     expand = c(0.01, 0.01),
                     limits = c(0, 25)) +
  ylab('Number of Gain or Loss Events (Thousands)') +
  xlab('phastCons Element')

ggsave('FigS21_gain_loss.pdf', gain_loss_plot, height = 1600, width = 2200, units = 'px')




