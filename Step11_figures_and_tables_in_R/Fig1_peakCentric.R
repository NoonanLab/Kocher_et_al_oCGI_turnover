# 10/19/22
# singleSpecies peakCentric plots for Fig 1
# https://clauswilke.com/blog/2016/06/13/reading-and-combining-many-tidy-data-files-in-r/

require(ggplot2)
require(tidyverse)
require(cowplot)
require(patchwork)
theme_set(theme_cowplot())

setwd('/Users/acadiak/Desktop/CGI/Figures/Mar_2023')

################################################################################
## WORK WITH PEAK-CENTRIC SUMMARY TABLES FOR MAKING FIG 1 C/D/E and FIG S7-12 ##
################################################################################

filePathToData <- '/Users/acadiak/Desktop/CGI/singleSpecies/peakCentric/summaryFiles_singleSpeciesPeakCentric/'

# store file names in a list
files <- dir(path = filePathToData, pattern = "*.txt")

# read in as nested dataframe
data <- tibble(filename = files) %>%
  mutate(file_contents = map(filename,
                             function(x) read_tsv(file.path(filePathToData, x)))
         )

# unnest to generate full table, but with file name as a column
unnested <- unnest(data, cols = c(file_contents))

# store levels for age categories
ageLevels <- c('Vertebrata', 'Gnathostomata', 'Tetrapoda', 'Amniota', 'Mammalia',
               'Theria', 'Eutheria', 'Primate', 'Ape', 'Human', 'None')

# modify variable names using various tidy commands
preFilter <- unnested %>% 
  # separate file names into new columns using separate
  separate(filename,
           c('species','tissue','mark','timePoint'),
           sep = '_',
           fill = 'right') %>%
  # change brain and limb to developing brain and developing limb in Noonan data (using the fact that only Noonan data has timePoint variable)
  mutate(newTissue = if_else(is.na(timePoint), tissue, paste('dev', tissue, sep=''))) %>%
  mutate(newTissue = factor(newTissue,
                            levels = c('brain', 'liver', 'muscle', 'testis', 'devbrain', 'devlimb'),
                            labels = c('Brain', 'Liver', 'Muscle', 'Testis', 'Dev Brain', 'Dev Limb'))) %>%
  # modify timePoint labels and order factors
  mutate(timePoint = factor(timePoint, 
                      levels = c('0.txt', '1.txt', '2.txt', '3.txt'), 
                      labels = c('0', '1', '2', '3'))) %>%
  # modify mark labels and order factors
  mutate(mark = factor(mark,
                       levels = c('H3K4me3.txt', 'H3K27ac.txt', 'H3K4me1.txt', 'ac', 'me2'),
                       labels = c('H3K4me3', 'H3K27ac', 'H3K4me1', 'H3K27ac', 'H3K4me2'))) %>%
  # modify species names and order factors
  mutate(species = factor(species,
                       levels = c('hg19', 'rheMac10', 'calJac4', 'mm39', 'rn7', 'susScr11', 'canFam6', 'felCat9', 'equCab3', 'rheMac2', 'mm9'),
                       labels = c('Human', 'Rhesus', 'Marmoset', 'Mouse', 'Rat', 'Pig', 'Dog', 'Cat', 'Horse', 'Rhesus', 'Mouse'))) %>%
  # transform certain columns into factors instead of numeric
  mutate_at(.vars = vars(species, tissue, mark,
                         CGI, peakLiftsToHuman, noFeaturesInHuman, 
                         Peak_OldestSegment, CGI_OldestSegment), factor) %>%
  # modify names for age segments and order factors
  mutate(Peak_OldestSegment = factor(Peak_OldestSegment,
                                     levels = rev(ageLevels), 
                                     labels = rev(c(rep('Older Than Amniota', 3), ageLevels[4:length(ageLevels)])))) %>%
  mutate(CGI_OldestSegment = factor(CGI_OldestSegment,
                                     levels = rev(ageLevels), 
                                     labels = rev(c(rep('Older Than Amniota', 3), ageLevels[4:length(ageLevels)]))))

# filter to only sites that lift to human and have no feature in human
table <- preFilter %>%
  filter(peakLiftsToHuman=='1'&
           noFeaturesInHuman=='1')

# write column newTissue to tissue
table$tissue <- table$newTissue
table <- table %>% dplyr::select(-newTissue)

# factor variables
table$CGI <- factor(table$CGI, levels = c(1,0))
table$mark <- factor(table$mark, levels = c('H3K4me3', 'H3K27ac', 'H3K4me1', 'H3K4me2'))
table$species <- factor(table$species, levels = c('Human', 'Rhesus', 'Marmoset', 'Mouse', 'Rat',
                                                  'Pig', 'Dog', 'Cat', 'Horse'))

# add label column combining mark and CGI
table$label <- paste(table$mark, table$CGI, sep = '_')
table$label <- factor(table$label, levels = c('H3K4me3_1', 'H3K4me3_0', 'H3K27ac_1', 'H3K27ac_0', 'H3K4me1_1', 'H3K4me1_0', 'H3K4me2_1', 'H3K4me2_0'))

# add length in thousands
table$length.in.thousands <- table$lengthPeak_Roller / 1000

# store color palette
cp.marks.alternateShades <- c('darkorange1', 'tan1', 
                              'palegreen4', 'darkseagreen3',
                              'mediumpurple4', 'mediumpurple2',
                              'darkgoldenrod1', 'gold')

# define function for getting 90% confidence interval, 25% quantile, median, 75% quantile
quantiles_90 <- function(x) {
  r <- quantile(x, probs=c(0.10, 0.25, 0.5, 0.75, 0.90))
  names(r) <- c("ymin", "lower", "middle", "upper", "ymax")
  r
}


################################################
# Peaks with oCGIs - Figure S7 (originally S4) #
################################################

peakCentric <- table %>% 
  group_by(species, tissue, mark, CGI) %>%
  dplyr::count() %>%
  mutate(tissue = factor(tissue,
                  levels = c('Brain', 'Liver', 'Muscle', 'Testis', 'Dev Brain', 'Dev Limb'),
                  labels = c('B', 'L', 'M', 'T', 'DB', 'DL'))) %>%
  pivot_wider(names_from = CGI, values_from = n) %>% tibble()
colnames(peakCentric) <- c('species', 'tissue', 'mark', 'n_CGI', 'n_noCGI', 'observed')
peakCentric$observed <- peakCentric$'n_CGI' / (peakCentric$'n_CGI' + peakCentric$'n_noCGI')
peakCentric$observed.pct <- peakCentric$observed * 100

# add extra entries to make plotting work
peakCentric <- peakCentric %>%
  add_row(species = 'Human', tissue = 'B', mark = 'H3K4me3',
          n_CGI = 0.0, n_noCGI = 0.0, observed = 0.0) %>%
  add_row(species = 'Human', tissue = 'B', mark = 'H3K4me1',
          n_CGI = 0.0, n_noCGI = 0.0, observed = 0.0)
for (species in c('Marmoset', 'Rat', 'Pig', 'Dog', 'Cat', 'Horse')) {
  peakCentric <- peakCentric %>%
    add_row(species = species, tissue = 'DB', mark = 'H3K4me2',
            n_CGI = 0.0, n_noCGI = 0.0, observed = 0.0)
}

# factor variables to get plots in the right order
peakCentric$species <- factor(peakCentric$species, levels = c('Human', 'Rhesus', 'Marmoset', 'Mouse', 'Rat', 'Pig', 'Dog', 'Cat', 'Horse'))
peakCentric$tissue <- factor(peakCentric$tissue, levels = c('B', 'L', 'M', 'T', 'DB', 'DL'))
peakCentric$mark <- factor(peakCentric$mark, levels = c('H3K4me3', 'H3K27ac', 'H3K4me1', 'H3K4me2'))

cp.marks <- c('darkorange1', 'palegreen4', 'mediumpurple4', 'darkgoldenrod1', 'gray30')

function_barplots_peakCentric <- function(table, inMark, color, ylimit) {
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
          strip.text = element_text(size = 16),
          panel.spacing.y = unit(1, 'lines')) +
    scale_y_continuous(labels = scales::percent,
                       limits = c(0, ylimit*1.1),
                       breaks = c(0, ylimit / 2, ylimit),
                       expand = c(0, 0)) +
    labs(y = '',
         x = '')
  return(x)
}

Peak_A <- function_barplots_peakCentric(peakCentric, 'H3K4me3', cp.marks[1], 0.6)
Peak_B <- function_barplots_peakCentric(peakCentric, 'H3K27ac', cp.marks[2], 0.3)
Peak_C <- function_barplots_peakCentric(peakCentric, 'H3K4me1', cp.marks[3], 0.3)
Peak_D <- function_barplots_peakCentric(peakCentric, 'H3K4me2', cp.marks[4], 0.3)

Peak_centric <- Peak_A + Peak_B + Peak_C + Peak_D + plot_layout(ncol = 4, widths = c(4, 6, 4, 1))

ggsave('FigS4_allSpecies_Peak-centric.pdf', Peak_centric, height = 2900, width = 2300, units = 'px')

# get exact percentages for the text
peakCentric %>%
  filter(species == 'Rhesus') %>%
  group_by(mark) %>%
  summarize(min = min(observed.pct),
            max = max(observed.pct))

# mark      min                 max
# H3K4me3   21.9 (testis)       49.5 (liver)
# H3K27ac   4.46 (muscle)       12.3 (dev brain)
# H3K4me1   4.96 (muscle)       7.82 (testis)
# H3K4me2   11.0 (dev brain)    11.0 (dev brain)


#####################################
### RPKM - Fig S8 (originally S5) ###
#####################################

# perform Wilcoxon test across every pair in the dataset
# (with CGI vs without CGI for all species x marks x tissues)
wilcox.RPKM <- table %>%
  group_by(species, mark, tissue) %>%
  do(w = wilcox.test(RPKM~CGI, data=., paired = F)) %>%
  summarize(species, mark, tissue, Wilcox = w$p.value)

# adjust p values
wilcox.RPKM$Wilcox.adj <- p.adjust(wilcox.RPKM$Wilcox, method = 'BH')

# factor variables
table$CGI <- factor(table$CGI, levels = c(1,0))
table$mark <- factor(table$mark, levels = c('H3K4me3', 'H3K27ac', 'H3K4me1', 'H3K4me2'))
table$species <- factor(table$species, levels = c('Human', 'Rhesus', 'Marmoset', 'Mouse', 'Rat',
                                                  'Pig', 'Dog', 'Cat', 'Horse'))

table$label <- paste(table$mark, table$CGI, sep = '_')
table$label <- factor(table$label, levels = c('H3K4me3_1', 'H3K4me3_0', 'H3K27ac_1', 'H3K27ac_0', 'H3K4me1_1', 'H3K4me1_0', 'H3K4me2_1', 'H3K4me2_0'))

cp.marks.alternateShades <- c('darkorange1', 'tan1', 
                              'palegreen4', 'darkseagreen3',
                              'mediumpurple4', 'mediumpurple2',
                              'darkgoldenrod1', 'gold')

# add info for significance lines to wilcox.RPKM
wilcox.RPKM$label <- paste0(wilcox.RPKM$mark, '_1')
wilcox.RPKM$x <- 0
wilcox.RPKM$x[wilcox.RPKM$tissue=='Brain'] <- 0.8
wilcox.RPKM$x[wilcox.RPKM$tissue=='Liver'] <- 1.8
wilcox.RPKM$x[wilcox.RPKM$tissue=='Muscle'] <- 2.8
wilcox.RPKM$x[wilcox.RPKM$tissue=='Testis'] <- 3.8
wilcox.RPKM$xend <- 0
wilcox.RPKM$xend[wilcox.RPKM$tissue=='Brain'] <- 1.2
wilcox.RPKM$xend[wilcox.RPKM$tissue=='Liver'] <- 2.2
wilcox.RPKM$xend[wilcox.RPKM$tissue=='Muscle'] <- 3.2
wilcox.RPKM$xend[wilcox.RPKM$tissue=='Testis'] <- 4.2
wilcox.RPKM$y <- 0
wilcox.RPKM$y[wilcox.RPKM$mark=='H3K4me3'] <- 22
wilcox.RPKM$y[wilcox.RPKM$mark=='H3K27ac'] <- 12
wilcox.RPKM$y[wilcox.RPKM$mark=='H3K4me1'] <- 12

wilcox.RPKM$x[wilcox.RPKM$tissue=='Dev Brain'] <- 0.8
wilcox.RPKM$x[wilcox.RPKM$tissue=='Dev Limb'] <- 1.8
wilcox.RPKM$xend[wilcox.RPKM$tissue=='Dev Brain'] <- 1.2
wilcox.RPKM$xend[wilcox.RPKM$tissue=='Dev Limb'] <- 2.2
wilcox.RPKM[(wilcox.RPKM$tissue=='Dev Brain'|wilcox.RPKM$tissue=='Dev Limb')&wilcox.RPKM$mark=='H3K27ac',]$y <- 22
wilcox.RPKM[(wilcox.RPKM$tissue=='Dev Brain'|wilcox.RPKM$tissue=='Dev Limb')&wilcox.RPKM$mark=='H3K4me2',]$y <- 8

wilcox.RPKM$yend <- wilcox.RPKM$y

# add info for significance stars to wilcox.RPKM
wilcox.RPKM$stars <- NA
wilcox.RPKM$stars[wilcox.RPKM$Wilcox.adj < 0.05] <- '*'
wilcox.RPKM$stars.x <- (wilcox.RPKM$x + wilcox.RPKM$xend) / 2

# drop rows without a star to prevent line from being plotted
wilcox.RPKM <- wilcox.RPKM[wilcox.RPKM$Wilcox.adj < 0.05,]

#stars <- tibble(species = c('Rhesus', 'Rhesus'),mark = c('H3K4me3', 'H3K4me3'),label = c('H3K4me3_1', 'H3K4me3_1'),x = c(1, 2),y = c(1.4, 1.4),stars = c('*', '*'))

# define function for getting 90% confidence interval, 25% quantile, median, 75% quantile
quantiles_90 <- function(x) {
  r <- quantile(x, probs=c(0.10, 0.25, 0.5, 0.75, 0.90))
  names(r) <- c("ymin", "lower", "middle", "upper", "ymax")
  r
}

# plot adult - Roller tissues, species Rhesus / Marmoset / Mouse / Rat
RPKM.adult.1 <- table %>%
  filter(tissue %in% c('Brain', 'Liver', 'Muscle', 'Testis'),
         species %in% c('Rhesus', 'Marmoset', 'Mouse', 'Rat')) %>%
  ggplot(aes(x = tissue, y = RPKM, fill = label)) +
  facet_grid(~factor(mark, levels = c('H3K4me3', 'H3K27ac', 'H3K4me1')) ~ 
               ~factor(species, levels = c('Rhesus', 'Marmoset', 'Mouse', 'Rat',
                                           'Pig', 'Dog', 'Cat', 'Horse')), 
             scales = 'free') +
  stat_summary(fun.data = quantiles_90, 
               geom='boxplot', 
               position=position_dodge2(preserve = 'single'), lwd = 0.25) +
  theme_bw() +
  scale_fill_manual(values = cp.marks.alternateShades[c(1:6)], guide = 'none') +
  theme(panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        strip.background = element_blank(),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 12),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        strip.text = element_text(size = 12)) +
  labs(y = 'Reads Per Kilobase Per\nMillion (RPKM) in Peak',
       x = '',
       size = 12) +
  expand_limits(y = 2) +
  geom_segment(data = wilcox.RPKM %>% filter(tissue %in% c('Brain', 'Liver', 'Muscle', 'Testis'),
                                               species %in% c('Rhesus', 'Marmoset', 'Mouse', 'Rat')), 
               aes(x = x, xend = xend, 
                   y = y, yend = yend), 
               inherit.aes = F,
               lwd = 0.5) +
  geom_text(data = wilcox.RPKM %>% filter(tissue %in% c('Brain', 'Liver', 'Muscle', 'Testis'),
                                            species %in% c('Rhesus', 'Marmoset', 'Mouse', 'Rat')),
            aes(x = stars.x , y = y + 0.1, label = stars), size = 6)

# print to output
ggsave('FigS5_partA1_RPKM_.pdf', RPKM.adult.1, height = 1100, width = 1800, units = 'px')

# same for Pig/Dog/Cat/Horse
RPKM.adult.2 <- table %>%
  filter(tissue %in% c('Brain', 'Liver', 'Muscle', 'Testis'),
         species %in% c('Pig', 'Dog', 'Cat', 'Horse')) %>%
  ggplot(aes(x = tissue, y = RPKM, fill = label)) +
  facet_grid(~factor(mark, levels = c('H3K4me3', 'H3K27ac', 'H3K4me1')) ~ 
               ~factor(species, levels = c('Rhesus', 'Marmoset', 'Mouse', 'Rat',
                                           'Pig', 'Dog', 'Cat', 'Horse')), 
             scales = 'free') +
  stat_summary(fun.data = quantiles_90, 
               geom='boxplot', 
               position=position_dodge2(preserve = 'single'), lwd = 0.25) +
  theme_bw() +
  scale_fill_manual(values = cp.marks.alternateShades[c(1:6)], guide = 'none') +
  theme(panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        strip.background = element_blank(),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 12),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        strip.text = element_text(size = 12)) +
  labs(y = 'Reads Per Kilobase Per\nMillion (RPKM) in Peak',
       x = '',
       size = 12) +
  expand_limits(y = 2) +
  geom_segment(data = wilcox.RPKM %>% filter(tissue %in% c('Brain', 'Liver', 'Muscle', 'Testis'),
                                             species %in% c('Pig', 'Dog', 'Cat', 'Horse')), 
               aes(x = x, xend = xend, 
                   y = y, yend = yend), 
               inherit.aes = F,
               lwd = 0.5) +
  geom_text(data = wilcox.RPKM %>% filter(tissue %in% c('Brain', 'Liver', 'Muscle', 'Testis'),
                                          species %in% c('Pig', 'Dog', 'Cat', 'Horse')),
            aes(x = stars.x , y = y + 0.1, label = stars), size = 6)

# print to output
ggsave('FigS5_partA2_RPKM_.pdf', RPKM.adult.2, height = 1100, width = 1800, units = 'px')

# same for Dev Brain / Dev limb H3K27ac
RPKM.dev.1 <- table %>%
  filter(tissue %in% c('Dev Brain', 'Dev Limb'),
         mark %in% c('H3K27ac')) %>%
  ggplot(aes(x = tissue, y = RPKM, fill = label)) +
  facet_grid(~factor(mark, levels = c('H3K27ac', 'H3K4me2')) ~ 
               ~factor(species, levels = c('Human', 'Rhesus', 'Mouse')), 
             scales = 'free') +
  stat_summary(fun.data = quantiles_90, 
               geom='boxplot', 
               position=position_dodge2(preserve = 'single'), lwd = 0.25) +
  theme_bw() +
  scale_fill_manual(values = cp.marks.alternateShades[c(3:4)], guide = 'none') +
  theme(panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        strip.background = element_blank(),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 12),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        strip.text = element_text(size = 12)) +
  labs(y = 'Average bigWig\nSignal Across Peak',
       x = '',
       size = 12) +
  expand_limits(y = 2) +
  geom_segment(data = wilcox.RPKM %>% filter(tissue %in% c('Dev Brain', 'Dev Limb'),
                                             mark %in% c('H3K27ac')), 
               aes(x = x, xend = xend, 
                   y = y, yend = yend), 
               inherit.aes = F,
               lwd = 0.5) +
  geom_text(data = wilcox.RPKM %>% filter(tissue %in% c('Dev Brain', 'Dev Limb'),
                                          mark %in% c('H3K27ac')),
            aes(x = stars.x , y = y + 0.1, label = stars), size = 6)

# print to output
ggsave('FigS5_partB1_RPKM_.pdf', RPKM.dev.1, height = 650, width = 800, units = 'px')

# same for Dev Brain H3K4me2
RPKM.dev.2 <- table %>%
  filter(tissue %in% c('Dev Brain'),
         mark %in% c('H3K4me2')) %>%
  ggplot(aes(x = tissue, y = RPKM, fill = label)) +
  facet_grid(~factor(mark, levels = c('H3K27ac', 'H3K4me2')) ~ 
               ~factor(species, levels = c('Human', 'Rhesus', 'Mouse')), 
             scales = 'free') +
  stat_summary(fun.data = quantiles_90, 
               geom='boxplot', 
               position=position_dodge2(preserve = 'single'), lwd = 0.25) +
  theme_bw() +
  scale_fill_manual(values = cp.marks.alternateShades[c(7:8)], guide = 'none') +
  theme(panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        strip.background = element_blank(),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 12),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        strip.text = element_text(size = 12)) +
  labs(y = 'Average bigWig\nSignal Across Peak',
       x = '',
       size = 12) +
  expand_limits(y = 2) +
  geom_segment(data = wilcox.RPKM %>% filter(tissue %in% c('Dev Brain'),
                                             mark %in% c('H3K4me2')), 
               aes(x = x, xend = xend, 
                   y = y, yend = yend), 
               inherit.aes = F,
               lwd = 0.5) +
  geom_text(data = wilcox.RPKM %>% filter(tissue %in% c('Dev Brain'),
                                          mark %in% c('H3K4me2')),
            aes(x = stars.x , y = y + 0.1, label = stars), size = 6)

# print to output
ggsave('FigS5_partB2_RPKM_.pdf', RPKM.dev.2, height = 650, width = 600, units = 'px')


################################################
###### PEAK LENGTH Fig S9 (originally S6) ######
################################################

# length of peaks with oCGIs vs peaks without oCGIs

# perform Wilcoxon test across every pair in the dataset
# (with CGI vs without CGI for all species x marks x tissues)
wilcox.length <- table %>%
  group_by(species, mark, tissue) %>%
  do(w = wilcox.test(lengthPeak_Roller~CGI, data=., paired = F)) %>%
  summarize(species, mark, tissue, Wilcox = w$p.value)

# adjust p values
wilcox.length$Wilcox.adj <- p.adjust(wilcox.length$Wilcox, method = 'BH')

# add info for significance lines to wilcox.length
wilcox.length$label <- paste0(wilcox.length$mark, '_1')
wilcox.length$x <- 0
wilcox.length$x[wilcox.length$tissue=='Brain'] <- 0.8
wilcox.length$x[wilcox.length$tissue=='Liver'] <- 1.8
wilcox.length$x[wilcox.length$tissue=='Muscle'] <- 2.8
wilcox.length$x[wilcox.length$tissue=='Testis'] <- 3.8
wilcox.length$xend <- 0
wilcox.length$xend[wilcox.length$tissue=='Brain'] <- 1.2
wilcox.length$xend[wilcox.length$tissue=='Liver'] <- 2.2
wilcox.length$xend[wilcox.length$tissue=='Muscle'] <- 3.2
wilcox.length$xend[wilcox.length$tissue=='Testis'] <- 4.2
wilcox.length$y <- 0
wilcox.length$y[wilcox.length$mark=='H3K4me3'] <- 1.3
wilcox.length$y[wilcox.length$mark=='H3K27ac'] <- 1.9
wilcox.length$y[wilcox.length$mark=='H3K4me1'] <- 2.5

wilcox.length$x[wilcox.length$tissue=='Dev Brain'] <- 0.8
wilcox.length$x[wilcox.length$tissue=='Dev Limb'] <- 1.8
wilcox.length$xend[wilcox.length$tissue=='Dev Brain'] <- 1.2
wilcox.length$xend[wilcox.length$tissue=='Dev Limb'] <- 2.2
wilcox.length[(wilcox.length$tissue=='Dev Brain'|wilcox.length$tissue=='Dev Limb')&wilcox.length$mark=='H3K27ac',]$y <- 4.5
wilcox.length[(wilcox.length$tissue=='Dev Brain'|wilcox.length$tissue=='Dev Limb')&wilcox.length$mark=='H3K4me2',]$y <- 3.1

# shift down for Pig/Dog/Cat/Horse H3K27ac and H3K4me1
wilcox.length[(wilcox.length$species=='Pig'|wilcox.length$species=='Dog'|wilcox.length$species=='Cat'|wilcox.length$species=='Horse')&wilcox.length$mark=='H3K27ac',]$y <- 1.4
wilcox.length[(wilcox.length$species=='Pig'|wilcox.length$species=='Dog'|wilcox.length$species=='Cat'|wilcox.length$species=='Horse')&wilcox.length$mark=='H3K4me1',]$y <- 1.9

wilcox.length$yend <- wilcox.length$y

# add info for significance stars to wilcox.length
wilcox.length$stars <- NA
wilcox.length$stars[wilcox.length$Wilcox.adj < 0.05] <- '*'
wilcox.length$stars.x <- (wilcox.length$x + wilcox.length$xend) / 2

# drop rows without a star to prevent line from being plotted
wilcox.length <- wilcox.length[wilcox.length$Wilcox.adj < 0.05,]

#### plots ####

# plot adult - Roller tissues, species Rhesus / Marmoset / Mouse / Rat
length.adult.1 <- table %>%
  filter(tissue %in% c('Brain', 'Liver', 'Muscle', 'Testis'),
         species %in% c('Rhesus', 'Marmoset', 'Mouse', 'Rat')) %>%
  ggplot(aes(x = tissue, y = length.in.thousands, fill = label)) +
  facet_grid(~factor(mark, levels = c('H3K4me3', 'H3K27ac', 'H3K4me1')) ~ 
               ~factor(species, levels = c('Rhesus', 'Marmoset', 'Mouse', 'Rat',
                                           'Pig', 'Dog', 'Cat', 'Horse')), 
             scales = 'free') +
  stat_summary(fun.data = quantiles_90, 
               geom='boxplot', 
               position=position_dodge2(preserve = 'single'), lwd = 0.25) +
  theme_bw() +
  scale_fill_manual(values = cp.marks.alternateShades[c(1:6)], guide = 'none') +
  theme(panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        strip.background = element_blank(),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 12),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        strip.text = element_text(size = 12)) +
  labs(y = 'Length of Peak \n(in Thousands of bp)',
       x = '',
       size = 12) +
  geom_segment(data = wilcox.length %>% filter(tissue %in% c('Brain', 'Liver', 'Muscle', 'Testis'),
                                               species %in% c('Rhesus', 'Marmoset', 'Mouse', 'Rat')), 
               aes(x = x, xend = xend, 
                   y = y, yend = yend), 
               inherit.aes = F,
               lwd = 0.5) +
  geom_text(data = wilcox.length %>% filter(tissue %in% c('Brain', 'Liver', 'Muscle', 'Testis'),
                                            species %in% c('Rhesus', 'Marmoset', 'Mouse', 'Rat')),
            aes(x = stars.x , y = y + 0.1, label = stars), size = 6)

# print to output
ggsave('FigS6_partA1_peakLength.pdf', length.adult.1, height = 1100, width = 1800, units = 'px')

# repeat for second 4 species (pig / dog / cat / horse)
length.adult.2 <- table %>%
  filter(tissue %in% c('Brain', 'Liver', 'Muscle', 'Testis'),
         species %in% c('Pig', 'Dog', 'Cat', 'Horse')) %>%
  ggplot(aes(x = tissue, y = length.in.thousands, fill = label)) +
  facet_grid(~factor(mark, levels = c('H3K4me3', 'H3K27ac', 'H3K4me1')) ~ 
               ~factor(species, levels = c('Rhesus', 'Marmoset', 'Mouse', 'Rat',
                                           'Pig', 'Dog', 'Cat', 'Horse')), 
             scales = 'free') +
  stat_summary(fun.data = quantiles_90, 
               geom='boxplot', 
               position=position_dodge2(preserve = 'single'), lwd = 0.25) +
  theme_bw() +
  scale_fill_manual(values = cp.marks.alternateShades[c(1:6)], guide = 'none') +
  theme(panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        strip.background = element_blank(),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 12),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        strip.text = element_text(size = 12)) +
  labs(y = 'Length of Peak \n(in Thousands of bp)',
       x = '',
       size = 12) +
  geom_segment(data = wilcox.length %>% filter(tissue %in% c('Brain', 'Liver', 'Muscle', 'Testis'),
                                               species %in% c('Pig', 'Dog', 'Cat', 'Horse')), 
               aes(x = x, xend = xend, 
                   y = y, yend = yend), 
               inherit.aes = F,
               lwd = 0.5) +
  geom_text(data = wilcox.length %>% filter(tissue %in% c('Brain', 'Liver', 'Muscle', 'Testis'),
                                            species %in% c('Pig', 'Dog', 'Cat', 'Horse')),
            aes(x = stars.x , y = y + 0.1, label = stars), size = 6)

# print to output
ggsave('FigS6_partA2_peakLength.pdf', length.adult.2, height = 1100, width = 1800, units = 'px')


# Noonan developing data
# first H3K27ac in both dev brain and dev limb
length.dev.1 <- table %>%
  filter(tissue %in% c('Dev Brain', 'Dev Limb'),
         mark %in% c('H3K27ac')) %>%
  ggplot(aes(x = tissue, y = length.in.thousands, fill = label)) +
  facet_grid(~factor(mark, levels = c('H3K27ac', 'H3K4me2')) ~ 
               ~factor(species, levels = c('Human', 'Rhesus', 'Mouse')), 
             scales = 'free') +
  stat_summary(fun.data = quantiles_90, 
               geom='boxplot', 
               position=position_dodge2(preserve = 'single'), lwd = 0.25) +
  theme_bw() +
  scale_fill_manual(values = cp.marks.alternateShades[c(3:4)], guide = 'none') +
  theme(panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        strip.background = element_blank(),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 12),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        strip.text = element_text(size = 12)) +
  labs(y = 'Length of Peak \n(in Thousands of bp)',
       x = '',
       size = 12) +
  geom_segment(data = wilcox.length %>% filter(tissue %in% c('Dev Brain', 'Dev Limb'),
                                               mark %in% c('H3K27ac')), 
               aes(x = x, xend = xend, 
                   y = y, yend = yend), 
               inherit.aes = F,
               lwd = 0.5) +
  geom_text(data = wilcox.length %>% filter(tissue %in% c('Dev Brain', 'Dev Limb'),
                                            mark %in% c('H3K27ac')),
            aes(x = stars.x , y = y + 0.1, label = stars), size = 6)

# print to output
ggsave('FigS6_partB1_peakLength_.pdf', length.dev.1, height = 650, width = 800, units = 'px')

# dev brain H3K4me2
length.dev.2 <- table %>%
  filter(tissue %in% c('Dev Brain'),
         mark %in% c('H3K4me2')) %>%
  ggplot(aes(x = tissue, y = length.in.thousands, fill = label)) +
  facet_grid(~factor(mark, levels = c('H3K27ac', 'H3K4me2')) ~ 
               ~factor(species, levels = c('Human', 'Rhesus', 'Mouse')), 
             scales = 'free') +
  stat_summary(fun.data = quantiles_90, 
               geom='boxplot', 
               position=position_dodge2(preserve = 'single'), lwd = 0.25) +
  theme_bw() +
  scale_fill_manual(values = cp.marks.alternateShades[c(7:8)], guide = 'none') +
  theme(panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        strip.background = element_blank(),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 12),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        strip.text = element_text(size = 12)) +
  labs(y = 'Length of Peak \n(in Thousands of bp)',
       x = '',
       size = 12) +
  geom_segment(data = wilcox.length %>% filter(tissue %in% c('Dev Brain'),
                                               mark %in% c('H3K4me2')), 
               aes(x = x, xend = xend, 
                   y = y, yend = yend), 
               inherit.aes = F,
               lwd = 0.5) +
  geom_text(data = wilcox.length %>% filter(tissue %in% c('Dev Brain'),
                                            mark %in% c('H3K4me2')),
            aes(x = stars.x , y = y + 0.1, label = stars), size = 6)

# print to output
ggsave('FigS6_partB2_peakLength_.pdf', length.dev.2, height = 650, width = 600, units = 'px')



###########################################
### phastCons - Fig S10 (originally S7) ###
###########################################

# perform Wilcoxon test across every pair in the dataset
# (with CGI vs without CGI for all species x marks x tissues)
wilcox.phastCons <- table %>%
  group_by(species, mark, tissue) %>%
  do(w = wilcox.test(Peak_phastConsMaxLOD~CGI, data=., paired = F)) %>%
  summarize(species, mark, tissue, Wilcox = w$p.value)

# adjust p values
wilcox.phastCons$Wilcox.adj <- p.adjust(wilcox.phastCons$Wilcox, method = 'BH')

# factor variables
table$CGI <- factor(table$CGI, levels = c(1,0))
table$mark <- factor(table$mark, levels = c('H3K4me3', 'H3K27ac', 'H3K4me1', 'H3K4me2'))
table$species <- factor(table$species, levels = c('Human', 'Rhesus', 'Marmoset', 'Mouse', 'Rat',
                                                  'Pig', 'Dog', 'Cat', 'Horse'))

table$label <- paste(table$mark, table$CGI, sep = '_')
table$label <- factor(table$label, levels = c('H3K4me3_1', 'H3K4me3_0', 'H3K27ac_1', 'H3K27ac_0', 'H3K4me1_1', 'H3K4me1_0', 'H3K4me2_1', 'H3K4me2_0'))

cp.marks.alternateShades <- c('darkorange1', 'tan1', 
                              'palegreen4', 'darkseagreen3',
                              'mediumpurple4', 'mediumpurple2',
                              'darkgoldenrod1', 'gold')

# add info for significance lines to wilcox.phastCons
wilcox.phastCons$label <- paste0(wilcox.phastCons$mark, '_1')
wilcox.phastCons$x <- 0
wilcox.phastCons$x[wilcox.phastCons$tissue=='Brain'] <- 0.8
wilcox.phastCons$x[wilcox.phastCons$tissue=='Liver'] <- 1.8
wilcox.phastCons$x[wilcox.phastCons$tissue=='Muscle'] <- 2.8
wilcox.phastCons$x[wilcox.phastCons$tissue=='Testis'] <- 3.8
wilcox.phastCons$xend <- 0
wilcox.phastCons$xend[wilcox.phastCons$tissue=='Brain'] <- 1.2
wilcox.phastCons$xend[wilcox.phastCons$tissue=='Liver'] <- 2.2
wilcox.phastCons$xend[wilcox.phastCons$tissue=='Muscle'] <- 3.2
wilcox.phastCons$xend[wilcox.phastCons$tissue=='Testis'] <- 4.2
wilcox.phastCons$y <- 0
wilcox.phastCons$y[wilcox.phastCons$mark=='H3K4me3'] <- 9
wilcox.phastCons$y[wilcox.phastCons$mark=='H3K27ac'] <- 8
wilcox.phastCons$y[wilcox.phastCons$mark=='H3K4me1'] <- 9

wilcox.phastCons$x[wilcox.phastCons$tissue=='Dev Brain'] <- 0.8
wilcox.phastCons$x[wilcox.phastCons$tissue=='Dev Limb'] <- 1.8
wilcox.phastCons$xend[wilcox.phastCons$tissue=='Dev Brain'] <- 1.2
wilcox.phastCons$xend[wilcox.phastCons$tissue=='Dev Limb'] <- 2.2
wilcox.phastCons[(wilcox.phastCons$tissue=='Dev Brain'|wilcox.phastCons$tissue=='Dev Limb')&wilcox.phastCons$mark=='H3K27ac',]$y <- 10
wilcox.phastCons[(wilcox.phastCons$tissue=='Dev Brain'|wilcox.phastCons$tissue=='Dev Limb')&wilcox.phastCons$mark=='H3K4me2',]$y <- 10

wilcox.phastCons$yend <- wilcox.phastCons$y

# add info for significance stars to wilcox.phastCons
wilcox.phastCons$stars <- NA
wilcox.phastCons$stars[wilcox.phastCons$Wilcox.adj < 0.05] <- '*'
wilcox.phastCons$stars.x <- (wilcox.phastCons$x + wilcox.phastCons$xend) / 2

# drop rows without a star to prevent line from being plotted
wilcox.phastCons <- wilcox.phastCons[wilcox.phastCons$Wilcox.adj <= 0.05,]

# define function for getting 90% confidence interval, 25% quantile, median, 75% quantile
quantiles_90 <- function(x) {
  r <- quantile(x, probs=c(0.10, 0.25, 0.5, 0.75, 0.90))
  names(r) <- c("ymin", "lower", "middle", "upper", "ymax")
  r
}

# get median and upper quartile values for the text
table$log2maxPhastCons <- log2(1 + table$Peak_phastConsMaxLOD)
withCGI <- table %>%
  filter(species == 'Rhesus', tissue == 'Brain', label == 'H3K4me3_1') %>%
  select(log2maxPhastCons)
quantile(withCGI$log2maxPhastCons, probs=c(0.10, 0.25, 0.5, 0.75, 0.90))
#      10%      25%      50%      75%      90% 
# 0.000000 0.000000 4.807355 7.025132 9.257152 
withoutCGI <- table %>%
  filter(species == 'Rhesus', tissue == 'Brain', label == 'H3K4me3_0') %>%
  select(log2maxPhastCons)
quantile(withoutCGI$log2maxPhastCons, probs=c(0.10, 0.25, 0.5, 0.75, 0.90))
#  10%      25%      50%      75%      90% 
# 0.000000 0.000000 3.807355 5.247928 7.294621 

# plot adult - Roller tissues, species Rhesus / Marmoset / Mouse / Rat
phastCons.adult.1 <- table %>%
  filter(tissue %in% c('Brain', 'Liver', 'Muscle', 'Testis'),
         species %in% c('Rhesus', 'Marmoset', 'Mouse', 'Rat')) %>%
  ggplot(aes(x = tissue, y = log2(1+Peak_phastConsMaxLOD), fill = label)) +
  facet_grid(~factor(mark, levels = c('H3K4me3', 'H3K27ac', 'H3K4me1')) ~ 
               ~factor(species, levels = c('Rhesus', 'Marmoset', 'Mouse', 'Rat',
                                           'Pig', 'Dog', 'Cat', 'Horse')), 
             scales = 'free') +
  stat_summary(fun.data = quantiles_90, 
               geom='boxplot', 
               position=position_dodge2(preserve = 'single'), lwd = 0.25) +
  theme_bw() +
  scale_fill_manual(values = cp.marks.alternateShades[c(1:6)], guide = 'none') +
  theme(panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        strip.background = element_blank(),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 12),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        strip.text = element_text(size = 12)) +
  labs(y = 'log2(1 + Max phastCons\nLOD Score in Peak)',
       x = '',
       size = 12) +
  geom_segment(data = wilcox.phastCons %>% filter(tissue %in% c('Brain', 'Liver', 'Muscle', 'Testis'),
                                             species %in% c('Rhesus', 'Marmoset', 'Mouse', 'Rat')), 
               aes(x = x, xend = xend, 
                   y = y, yend = yend), 
               inherit.aes = F,
               lwd = 0.5) +
  geom_text(data = wilcox.phastCons %>% filter(tissue %in% c('Brain', 'Liver', 'Muscle', 'Testis'),
                                          species %in% c('Rhesus', 'Marmoset', 'Mouse', 'Rat')),
            aes(x = stars.x , y = y + 0.1, label = stars), size = 6)

# print to output
ggsave('FigS7_partA1_phastCons_.pdf', phastCons.adult.1, height = 1100, width = 1800, units = 'px')

# same for Pig/Dog/Cat/Horse
phastCons.adult.2 <- table %>%
  filter(tissue %in% c('Brain', 'Liver', 'Muscle', 'Testis'),
         species %in% c('Pig', 'Dog', 'Cat', 'Horse')) %>%
  ggplot(aes(x = tissue, y = log2(1+Peak_phastConsMaxLOD), fill = label)) +
  facet_grid(~factor(mark, levels = c('H3K4me3', 'H3K27ac', 'H3K4me1')) ~ 
               ~factor(species, levels = c('Rhesus', 'Marmoset', 'Mouse', 'Rat',
                                           'Pig', 'Dog', 'Cat', 'Horse')), 
             scales = 'free') +
  stat_summary(fun.data = quantiles_90, 
               geom='boxplot', 
               position=position_dodge2(preserve = 'single'), lwd = 0.25) +
  theme_bw() +
  scale_fill_manual(values = cp.marks.alternateShades[c(1:6)], guide = 'none') +
  theme(panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        strip.background = element_blank(),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 12),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        strip.text = element_text(size = 12)) +
  labs(y = 'log2(1 + Max phastCons\nLOD Score in Peak)',
       x = '',
       size = 12) +
  expand_limits(y = 2) +
  coord_cartesian(ylim = c(0, 10)) +
  geom_segment(data = wilcox.phastCons %>% filter(tissue %in% c('Brain', 'Liver', 'Muscle', 'Testis'),
                                                  species %in% c('Pig', 'Dog', 'Cat', 'Horse')), 
               aes(x = x, xend = xend, 
                   y = y, yend = yend), 
               inherit.aes = F,
               lwd = 0.5) +
  geom_text(data = wilcox.phastCons %>% filter(tissue %in% c('Brain', 'Liver', 'Muscle', 'Testis'),
                                               species %in% c('Pig', 'Dog', 'Cat', 'Horse')),
            aes(x = stars.x , y = y + 0.1, label = stars), size = 6)

# print to output
ggsave('FigS7_partA2_phastCons_.pdf', phastCons.adult.2, height = 1100, width = 1800, units = 'px')

# same for Dev Brain / Dev limb H3K27ac
phastCons.dev.1 <- table %>%
  filter(tissue %in% c('Dev Brain', 'Dev Limb'),
         mark %in% c('H3K27ac')) %>%
  ggplot(aes(x = tissue, y = log2(1+Peak_phastConsMaxLOD), fill = label)) +
  facet_grid(~factor(mark, levels = c('H3K27ac', 'H3K4me2')) ~ 
               ~factor(species, levels = c('Human', 'Rhesus', 'Mouse')), 
             scales = 'free') +
  stat_summary(fun.data = quantiles_90, 
               geom='boxplot', 
               position=position_dodge2(preserve = 'single'), lwd = 0.25) +
  theme_bw() +
  scale_fill_manual(values = cp.marks.alternateShades[c(3:4)], guide = 'none') +
  theme(panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        strip.background = element_blank(),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 12),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        strip.text = element_text(size = 12)) +
  labs(y = 'log2(1 + Max phastCons\nLOD Score in Peak)',
       x = '',
       size = 12) +
  expand_limits(y = 2) +
  geom_segment(data = wilcox.phastCons %>% filter(tissue %in% c('Dev Brain', 'Dev Limb'),
                                             mark %in% c('H3K27ac')), 
               aes(x = x, xend = xend, 
                   y = y, yend = yend), 
               inherit.aes = F,
               lwd = 0.5) +
  geom_text(data = wilcox.phastCons %>% filter(tissue %in% c('Dev Brain', 'Dev Limb'),
                                          mark %in% c('H3K27ac')),
            aes(x = stars.x , y = y + 0.1, label = stars), size = 6)

# print to output
ggsave('FigS7_partB1_phastCons_.pdf', phastCons.dev.1, height = 650, width = 800, units = 'px')

# same for Dev Brain H3K4me2
phastCons.dev.2 <- table %>%
  filter(tissue %in% c('Dev Brain'),
         mark %in% c('H3K4me2')) %>%
  ggplot(aes(x = tissue, y = log2(1+Peak_phastConsMaxLOD), fill = label)) +
  facet_grid(~factor(mark, levels = c('H3K27ac', 'H3K4me2')) ~ 
               ~factor(species, levels = c('Human', 'Rhesus', 'Mouse')), 
             scales = 'free') +
  stat_summary(fun.data = quantiles_90, 
               geom='boxplot', 
               position=position_dodge2(preserve = 'single'), lwd = 0.25) +
  theme_bw() +
  scale_fill_manual(values = cp.marks.alternateShades[c(7:8)], guide = 'none') +
  theme(panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        strip.background = element_blank(),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 12),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        strip.text = element_text(size = 12)) +
  labs(y = 'log2(1 + Max phastCons\nLOD Score in Peak)',
       x = '',
       size = 12) +
  expand_limits(y = 2) +
  geom_segment(data = wilcox.phastCons %>% filter(tissue %in% c('Dev Brain'),
                                             mark %in% c('H3K4me2')), 
               aes(x = x, xend = xend, 
                   y = y, yend = yend), 
               inherit.aes = F,
               lwd = 0.5) +
  geom_text(data = wilcox.phastCons %>% filter(tissue %in% c('Dev Brain'),
                                          mark %in% c('H3K4me2')),
            aes(x = stars.x , y = y + 0.1, label = stars), size = 6)

# print to output
ggsave('FigS7_partB2_phastCons_.pdf', phastCons.dev.2, height = 650, width = 600, units = 'px')


#############################################
### sequence age - Fig 12 (originally S9) ###
#############################################

# For non-ape species: convert all ape/human age sequences to None
# For non-primate species: convert all primate/ape/human age sequences to None

table[table$species!='Human'&table$Peak_OldestSegment=='Human',]$Peak_OldestSegment <- 'None'
table[table$species!='Human'&table$Peak_OldestSegment=='Ape',]$Peak_OldestSegment <- 'None'

table[table$species!='Human'&table$species!='Rhesus'&table$species!='Marmoset'&table$Peak_OldestSegment=='Human',]$Peak_OldestSegment <- 'None'
table[table$species!='Human'&table$species!='Rhesus'&table$species!='Marmoset'&table$Peak_OldestSegment=='Ape',]$Peak_OldestSegment <- 'None'
table[table$species!='Human'&table$species!='Rhesus'&table$species!='Marmoset'&table$Peak_OldestSegment=='Primate',]$Peak_OldestSegment <- 'None'

# count categories for ADULT tissues (Roller)
ageCounts <- table %>%
  filter(tissue %in% c('Brain', 'Liver', 'Muscle', 'Testis')) %>%
  group_by(species, mark, tissue, CGI, Peak_OldestSegment) %>%
  count()

age.rhesus <- ageCounts %>%
    filter(species == 'Rhesus')  %>%
    ggplot(aes(x = CGI, y = n, fill = Peak_OldestSegment)) +
    facet_grid(mark ~ tissue) +
    geom_bar(position = 'fill', stat = 'identity') +
    theme_bw() +
    theme(panel.grid.minor.x = element_blank(),
          panel.grid.major.x = element_blank(),
          panel.grid.minor.y = element_blank(),
          strip.background = element_blank(),
          axis.text = element_text(size = 12),
          axis.title = element_text(size = 12),
          axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
          strip.text = element_text(size = 12),
          strip.text.y = element_blank(),
          legend.position = "none") +
    scale_y_continuous(labels = scales::percent,
                       breaks = c(0, 0.5, 1.0)) +
    scale_fill_manual(values = c('gray90', 'gray70', 'gray50', '#c994c7', '#df65b0', '#dd1c77','#980043')) +
    labs(y = 'Percent of Peaks with\nOldest Sequence of Each Age',
         x = '')
 
age.marm <- ageCounts %>%
  filter(species == 'Marmoset')  %>%
  ggplot(aes(x = CGI, y = n, fill = Peak_OldestSegment)) +
  facet_grid(mark ~ tissue) +
  geom_bar(position = 'fill', stat = 'identity') +
  theme_bw() +
  theme(panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        strip.background = element_blank(),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 12),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        strip.text.y = element_blank(),
        strip.text = element_text(size = 12),
        legend.position = "none") +
  scale_fill_manual(values = c('gray90', 'gray70', 'gray50', '#c994c7', '#df65b0', '#dd1c77','#980043')) +
  labs(y = '',
       x = '')

age.mouse <- ageCounts %>%
  filter(species == 'Mouse')  %>%
  ggplot(aes(x = CGI, y = n, fill = Peak_OldestSegment)) +
  facet_grid(mark ~ tissue) +
  geom_bar(position = 'fill', stat = 'identity') +
  theme_bw() +
  theme(panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        strip.background = element_blank(),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 12),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        strip.text.y = element_blank(),
        strip.text = element_text(size = 12),
        legend.position = "none") +
  scale_fill_manual(values = c('gray90', 'gray50', '#c994c7', '#df65b0', '#dd1c77','#980043')) +
  labs(y = '',
       x = '')

age.rat <- ageCounts %>%
  filter(species == 'Rat')  %>%
  ggplot(aes(x = CGI, y = n, fill = Peak_OldestSegment)) +
  facet_grid(mark ~ tissue) +
  geom_bar(position = 'fill', stat = 'identity') +
  theme_bw() +
  theme(panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        strip.background = element_blank(),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 12),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        strip.text = element_text(size = 12),
        legend.position = "none") +
  scale_fill_manual(values = c('gray90', 'gray50', '#c994c7', '#df65b0', '#dd1c77','#980043')) +
  labs(y = '',
       x = '')

age.pig <- ageCounts %>%
  filter(species == 'Pig')  %>%
  ggplot(aes(x = CGI, y = n, fill = Peak_OldestSegment)) +
  facet_grid(mark ~ tissue) +
  geom_bar(position = 'fill', stat = 'identity') +
  theme_bw() +
  theme(panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        strip.background = element_blank(),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 12),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        strip.text = element_text(size = 12),
        strip.text.y = element_blank(),
        legend.position = "none") +
  scale_y_continuous(labels = scales::percent,
                     breaks = c(0, 0.5, 1.0)) +
  scale_fill_manual(values = c('gray90', 'gray50', '#c994c7', '#df65b0', '#dd1c77','#980043')) +
  labs(y = 'Percent of Peaks with\nOldest Sequence of Each Age',
       x = '')

age.dog <- ageCounts %>%
  filter(species == 'Dog')  %>%
  ggplot(aes(x = CGI, y = n, fill = Peak_OldestSegment)) +
  facet_grid(mark ~ tissue) +
  geom_bar(position = 'fill', stat = 'identity') +
  theme_bw() +
  theme(panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        strip.background = element_blank(),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 12),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        strip.text.y = element_blank(),
        strip.text = element_text(size = 12),
        legend.position = "none") +
  scale_fill_manual(values = c('gray90', 'gray50', '#c994c7', '#df65b0', '#dd1c77','#980043')) +
  labs(y = '',
       x = '')

age.cat <- ageCounts %>%
  filter(species == 'Cat')  %>%
  ggplot(aes(x = CGI, y = n, fill = Peak_OldestSegment)) +
  facet_grid(mark ~ tissue) +
  geom_bar(position = 'fill', stat = 'identity') +
  theme_bw() +
  theme(panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        strip.background = element_blank(),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 12),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        strip.text.y = element_blank(),
        strip.text = element_text(size = 12),
        legend.position = "none") +
  scale_fill_manual(values = c('gray90', 'gray50', '#c994c7', '#df65b0', '#dd1c77','#980043')) +
  labs(y = '',
       x = '')

age.horse <- ageCounts %>%
  filter(species == 'Horse')  %>%
  ggplot(aes(x = CGI, y = n, fill = Peak_OldestSegment)) +
  facet_grid(mark ~ tissue) +
  geom_bar(position = 'fill', stat = 'identity') +
  theme_bw() +
  theme(panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        strip.background = element_blank(),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 12),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        strip.text = element_text(size = 12),
        legend.position = "none") +
  scale_fill_manual(values = c('gray90', 'gray50', '#c994c7', '#df65b0', '#dd1c77','#980043')) +
  labs(y = '',
       x = '')

age.adult <- (age.rhesus + age.marm + age.mouse + age.rat + plot_layout(ncol = 4)) / (age.pig + age.dog + age.cat + age.horse + plot_layout(ncol = 4))
ggsave('FigS9_partA_peakAge.pdf', age.adult, height = 1900, width = 2200, units = 'px')

###### COLOR KEY
# c('gray90', 'gray85', 'gray80', 'gray70', 'gray50', '#c994c7', '#df65b0', '#dd1c77','#980043')
# c('None',    'Human',  'Ape',   'Primate',  'Eutheria', 'Theria', 'Mammalia', 'Amniota', 'Older Than Amniota)

############
# count categories for DEVELOPING tissues (Noonan)
ageCounts <- table %>%
  filter(tissue %in% c('Dev Brain', 'Dev Limb')) %>%
  group_by(species, mark, tissue, CGI, Peak_OldestSegment) %>%
  count()

# developing H3K27ac

age.human.devK27ac <- ageCounts %>%
  filter(species == 'Human', mark == 'H3K27ac')  %>%
  ggplot(aes(x = CGI, y = n, fill = Peak_OldestSegment)) +
  facet_grid(mark ~ tissue) +
  geom_bar(position = 'fill', stat = 'identity') +
  theme_bw() +
  theme(panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        strip.background = element_blank(),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 12),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        strip.text = element_text(size = 12),
        strip.text.y = element_blank(),
        legend.position = "none") +
  scale_y_continuous(labels = scales::percent,
                     breaks = c(0, 0.5, 1.0)) +
  scale_fill_manual(values = c('gray90', 'gray83', 'gray77', 'gray70', 'gray50', '#c994c7', '#df65b0', '#dd1c77','#980043')) +
  labs(y = 'Percent of Peaks with\nOldest Sequence of Each Age',
       x = '')

age.rhesus.devK27ac <- ageCounts %>%
  filter(species == 'Rhesus', mark == 'H3K27ac')  %>%
  ggplot(aes(x = CGI, y = n, fill = Peak_OldestSegment)) +
  facet_grid(mark ~ tissue) +
  geom_bar(position = 'fill', stat = 'identity') +
  theme_bw() +
  theme(panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        strip.background = element_blank(),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 12),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        strip.text.y = element_blank(),
        strip.text = element_text(size = 12),
        legend.position = "none") +
  scale_fill_manual(values = c('gray90', 'gray70', 'gray50', '#c994c7', '#df65b0', '#dd1c77','#980043')) +
  labs(y = '',
       x = '')

age.mouse.devK27ac <- ageCounts %>%
  filter(species == 'Mouse', mark == 'H3K27ac')  %>%
  ggplot(aes(x = CGI, y = n, fill = Peak_OldestSegment)) +
  facet_grid(mark ~ tissue) +
  geom_bar(position = 'fill', stat = 'identity') +
  theme_bw() +
  theme(panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        strip.background = element_blank(),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 12),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        strip.text.y = element_blank(),
        strip.text = element_text(size = 12),
        legend.position = "none") +
  scale_fill_manual(values = c('gray90', 'gray50', '#c994c7', '#df65b0', '#dd1c77','#980043')) +
  labs(y = '',
       x = '')

###### COLOR KEY
# c('gray90', 'gray85', 'gray80', 'gray70', 'gray50', '#c994c7', '#df65b0', '#dd1c77','#980043')
# c('None',    'Human',  'Ape',   'Primate',  'Eutheria', 'Theria', 'Mammalia', 'Amniota', 'Older Than Amniota)

age.dev.H3K27ac <- age.human.devK27ac + age.rhesus.devK27ac + age.mouse.devK27ac + plot_layout(ncol = 3)
ggsave('FigS9_partB_peakAge.pdf', age.dev.H3K27ac, height = 550, width = 1000, units = 'px')

## developing H3K4me2

age.human.devK4me2 <- ageCounts %>%
  filter(species == 'Human', mark == 'H3K4me2')  %>%
  ggplot(aes(x = CGI, y = n, fill = Peak_OldestSegment)) +
  facet_grid(mark ~ tissue) +
  geom_bar(position = 'fill', stat = 'identity') +
  theme_bw() +
  theme(panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        strip.background = element_blank(),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 12),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        strip.text = element_text(size = 12),
        strip.text.y = element_blank(),
        legend.position = "none") +
  scale_y_continuous(labels = scales::percent,
                     breaks = c(0, 0.5, 1.0)) +
  scale_fill_manual(values = c('gray90', 'gray83', 'gray77', 'gray70', 'gray50', '#c994c7', '#df65b0', '#dd1c77','#980043')) +
  labs(y = 'Percent of Peaks with\nOldest Sequence of Each Age',
       x = '')

age.rhesus.devK4me2 <- ageCounts %>%
  filter(species == 'Rhesus', mark == 'H3K4me2')  %>%
  ggplot(aes(x = CGI, y = n, fill = Peak_OldestSegment)) +
  facet_grid(mark ~ tissue) +
  geom_bar(position = 'fill', stat = 'identity') +
  theme_bw() +
  theme(panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        strip.background = element_blank(),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 12),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        strip.text.y = element_blank(),
        strip.text = element_text(size = 12),
        legend.position = "none") +
  scale_fill_manual(values = c('gray90', 'gray70', 'gray50', '#c994c7', '#df65b0', '#dd1c77','#980043')) +
  labs(y = '',
       x = '')

age.mouse.devK4me2 <- ageCounts %>%
  filter(species == 'Mouse', mark == 'H3K4me2')  %>%
  ggplot(aes(x = CGI, y = n, fill = Peak_OldestSegment)) +
  facet_grid(mark ~ tissue) +
  geom_bar(position = 'fill', stat = 'identity') +
  theme_bw() +
  theme(panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        strip.background = element_blank(),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 12),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        strip.text.y = element_blank(),
        strip.text = element_text(size = 12),
        legend.position = "none") +
  scale_fill_manual(values = c('gray90', 'gray50', '#c994c7', '#df65b0', '#dd1c77','#980043')) +
  labs(y = '',
       x = '')

###### COLOR KEY
# c('gray90', 'gray85', 'gray80', 'gray70', 'gray50', '#c994c7', '#df65b0', '#dd1c77','#980043')
# c('None',    'Human',  'Ape',   'Primate',  'Eutheria', 'Theria', 'Mammalia', 'Amniota', 'Older Than Amniota)

age.dev.H3K4me2 <- age.human.devK4me2 + age.rhesus.devK4me2 + age.mouse.devK4me2 + plot_layout(ncol = 3)
ggsave('FigS9_partC_peakAge.pdf', age.dev.H3K4me2, height = 550, width = 800, units = 'px')


# plot single human dev plot to get legend

age.legend <- ageCounts %>%
  filter(species == 'Human', mark == 'H3K27ac')  %>%
  ggplot(aes(x = CGI, y = n, fill = Peak_OldestSegment)) +
  facet_grid(mark ~ tissue) +
  geom_bar(position = 'fill', stat = 'identity') +
  theme_bw() +
  theme(panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        strip.background = element_blank(),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 12),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        strip.text = element_text(size = 12),
        strip.text.y = element_blank()) +
  scale_y_continuous(labels = scales::percent,
                     breaks = c(0, 0.5, 1.0)) +
  scale_fill_manual(values = c('gray90', 'gray83', 'gray77', 'gray70', 'gray50', '#c994c7', '#df65b0', '#dd1c77','#980043')) +
  labs(y = 'Percent of Peaks with\nOldest Sequence of Each Age',
       x = '')
ggsave('FigS9_legend_peakAge.pdf', age.legend, height = 800, width = 1200, units = 'px')


###########################################################################
#### FIG 1 C/D/E and Fig S11 (originally S8) - rhesus brain as example ####
###########################################################################

# RPKM for Fig 1C
RPKM.signif <- data.frame(x=c(0.8, 1.8, 2.8),
                            xend=c(1.2, 2.2, 3.2),
                            y=rep(11, 3),
                            yend=rep(11,3),
                            stars=c('*','*','*'),
                            stars.x=c(1,2,3),
                            label=c('H3K4me3_1','H3K27ac_1','H3K4me1_1'))

RPKM.rhesus <- table %>%
  filter(tissue %in% c('Brain'),
         species %in% c('Rhesus')) %>%
  ggplot(aes(x = mark, y = RPKM, fill = label)) +
  stat_summary(fun.data = quantiles_90, 
               geom='boxplot', 
               position=position_dodge2(preserve = 'single'), lwd = 0.25) +
  theme_bw() +
  scale_fill_manual(values = cp.marks.alternateShades[c(1:6)], guide = 'none') +
  theme(panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        strip.background = element_blank(),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 12),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        strip.text = element_text(size = 12)) +
  labs(y = 'Reads Per Kilobase Per\nMillion (RPKM) in Peak',
       x = '',
       size = 12) +
  expand_limits(y = 2) +
  geom_segment(data = RPKM.signif, 
               aes(x = x, xend = xend, 
                   y = y, yend = yend),
               inherit.aes = F,
               lwd = 0.5) +
  geom_text(data = RPKM.signif,
            aes(x = stars.x , y = y + 0.1, label = stars), size = 6)
ggsave('Fig1C_peakRPKM_rhesus.pdf', RPKM.rhesus, height = 850, width = 700, units = 'px')

# phastCons Max LOD for Fig 1D
phastCons.maxLOD.signif <- data.frame(x=c(0.8, 1.8, 2.8),
                          xend=c(1.2, 2.2, 3.2),
                          y=rep(8, 3),
                          yend=rep(8,3),
                          stars=c('*','*','*'),
                          stars.x=c(1,2,3),
                          label=c('H3K4me3_1','H3K27ac_1','H3K4me1_1'))

phastCons.maxLOD.rhesus <- table %>%
  filter(tissue %in% c('Brain'),
         species %in% c('Rhesus')) %>%
  ggplot(aes(x = mark, y = log2(1+Peak_phastConsMaxLOD), fill = label)) +
  stat_summary(fun.data = quantiles_90, 
               geom='boxplot', 
               position=position_dodge2(preserve = 'single'), lwd = 0.25) +
  theme_bw() +
  scale_fill_manual(values = cp.marks.alternateShades[c(1:6)], guide = 'none') +
  theme(panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        strip.background = element_blank(),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 12),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        strip.text = element_text(size = 12)) +
  labs(y = 'log2(1 + Max phastCons\nLOD Score in Peak)',
       x = '',
       size = 12) +
  geom_segment(data = phastCons.maxLOD.signif, 
               aes(x = x, xend = xend, 
                   y = y, yend = yend), 
               inherit.aes = F,
               lwd = 0.5) +
  geom_text(data = phastCons.maxLOD.signif,
            aes(x = stars.x , y = y + 0.1, label = stars), size = 6)
ggsave('Fig1D_phastCons_maxLOD_rhesus.pdf', phastCons.maxLOD.rhesus, height = 850, width = 700, units = 'px')

# age distribution for Fig 1E

ageCounts <- table %>%
  filter(species %in% c('Rhesus'),
         tissue %in% c('Brain')) %>%
  group_by(species, mark, tissue, CGI, Peak_OldestSegment) %>%
  count()

age.rhesus.main <- ageCounts %>%
  filter(species == 'Rhesus')  %>%
  ggplot(aes(x = CGI, y = n, fill = Peak_OldestSegment)) +
  facet_grid(. ~ mark) +
  geom_bar(position = 'fill', stat = 'identity') +
  theme_bw() +
  theme(panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        strip.background = element_blank(),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 12),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        strip.text = element_text(size = 12),
        strip.text.y = element_blank()) +
  scale_y_continuous(labels = scales::percent,
                     breaks = c(0, 0.5, 1.0)) +
  scale_fill_manual(values = c('gray90', 'gray70', 'gray50', '#c994c7', '#df65b0', '#dd1c77','#980043')) +
  labs(y = 'Percent of Peaks with\nOldest Sequence of Each Age',
       x = '')

ggsave('Fig1E_age_rhesus.pdf', age.rhesus.main, height = 800, width = 1300, units = 'px')


# ALT PHASTCONS MEASURES FOR FIG S8
# DO WILCOXON TEST ON THE TWO BELOW, AND ADJUST BETWEEN THE TWO OF THEM

# do wilcoxon test for sum of LOD scores and percent Bp in phastCons element
wilcox.phastCons.sumLOD <- table %>%
  filter(tissue %in% c('Brain'),
         species %in% c('Rhesus')) %>%
  group_by(mark) %>%
  do(w = wilcox.test(Peak_phastConsSumLOD~CGI, data=., paired = F)) %>%
  summarize(mark, Wilcox = w$p.value)
wilcox.phastCons.sumLOD$test <- 'sumLOD'

table$Peak_phastConsPercentBp <- table$Peak_phastConsBp / table$lengthPeak_Human
wilcox.phastCons.percentBp <- table %>%
  filter(tissue %in% c('Brain'),
         species %in% c('Rhesus')) %>%
  group_by(mark) %>%
  do(w = wilcox.test(Peak_phastConsPercentBp~CGI, data=., paired = F)) %>%
  summarize(mark, Wilcox = w$p.value)
wilcox.phastCons.percentBp$test <- 'percentBp'

# paste together two tables and calculate adjusted p value
wilcox.alternatePhastCons <- full_join(wilcox.phastCons.sumLOD, wilcox.phastCons.percentBp)
wilcox.alternatePhastCons$Wilcox.adj <- p.adjust(wilcox.alternatePhastCons$Wilcox, method = 'BH')

# add columns with info for adding stars
wilcox.alternatePhastCons$x <- rep(c(0.8, 1.8, 2.8), 2)
wilcox.alternatePhastCons$xend <- wilcox.alternatePhastCons$x + 0.4
wilcox.alternatePhastCons$y <- c(9,9,9,0.35,0.35,0.35)
wilcox.alternatePhastCons$yend <- wilcox.alternatePhastCons$y
wilcox.alternatePhastCons$stars <- rep('*', 6)
wilcox.alternatePhastCons$stars.x <- wilcox.alternatePhastCons$x + 0.2
wilcox.alternatePhastCons$label <- paste0(wilcox.alternatePhastCons$mark, '_1')

# make plots
phastCons.sumLOD.rhesus <- table %>%
  filter(tissue %in% c('Brain'),
         species %in% c('Rhesus')) %>%
  ggplot(aes(x = mark, y = log2(1+Peak_phastConsSumLOD), fill = label)) +
  stat_summary(fun.data = quantiles_90, 
               geom='boxplot', 
               position=position_dodge2(preserve = 'single'), lwd = 0.25) +
  theme_bw() +
  scale_fill_manual(values = cp.marks.alternateShades[c(1:6)], guide = 'none') +
  theme(panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        strip.background = element_blank(),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 12),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        strip.text = element_text(size = 12)) +
  labs(y = 'log2(1 + Sum of phastCons\nLOD Scores in Peak)',
       x = '',
       size = 12) +
  geom_segment(data = wilcox.alternatePhastCons %>% filter(test == 'sumLOD'), 
               aes(x = x, xend = xend, 
                   y = y, yend = yend), 
               inherit.aes = F,
               lwd = 0.5) +
  geom_text(data = wilcox.alternatePhastCons %>% filter(test == 'sumLOD'),
            aes(x = stars.x , y = y + 0.1, label = stars), size = 6)
ggsave('FigS8_phastCons_sumLOD_rhesus.pdf', phastCons.sumLOD.rhesus, height = 850, width = 700, units = 'px')

# phastCons % Bp in phastCons element
phastCons.bpPercent.rhesus <- table %>%
  filter(tissue %in% c('Brain'),
         species %in% c('Rhesus')) %>%
  ggplot(aes(x = mark, y = Peak_phastConsPercentBp, fill = label)) +
  stat_summary(fun.data = quantiles_90, 
               geom='boxplot', 
               position=position_dodge2(preserve = 'single'), lwd = 0.25) +
  theme_bw() +
  scale_fill_manual(values = cp.marks.alternateShades[c(1:6)], guide = 'none') +
  theme(panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        strip.background = element_blank(),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 12),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        strip.text = element_text(size = 12)) +
  labs(y = 'Percent of Peak Overlapping\na phastCons Element',
       x = '',
       size = 12) +
  scale_y_continuous(labels = scales::percent) +
  geom_segment(data = wilcox.alternatePhastCons %>% filter(test == 'percentBp'), 
               aes(x = x, xend = xend, 
                   y = y, yend = yend), 
               inherit.aes = F,
               lwd = 0.5) +
  geom_text(data = wilcox.alternatePhastCons %>% filter(test == 'percentBp'),
            aes(x = stars.x , y = y + 0.01, label = stars), size = 6)
ggsave('FigS8_phastCons_bpPercent_rhesus.pdf', phastCons.bpPercent.rhesus, height = 850, width = 700, units = 'px')


