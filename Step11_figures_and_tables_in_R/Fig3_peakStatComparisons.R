# 11/10/22
# Purpose: make plots in Fig 3E, several supplements (peak conservation/age)
# for rhesus vs mouse H3K27ac peaks

require(ggplot2)
require(cowplot)
theme_set(theme_cowplot())
require(tidyverse)
require(patchwork)

setwd('/Users/acadiak/Desktop/CGI/Figures/Mar_2023')

# READ IN DATA FROM ALL SPECIES PAIRS

# store file names in a list
filePathToData <- '/Users/acadiak/Desktop/CGI/speciesPairs/CGIcentric/Roller_summaryFiles/'
files <- dir(path = filePathToData, pattern = "*.txt")

# read in as nested dataframe
dataRoller <- tibble(filename = files) %>%
  mutate(file_contents = map(filename,
                             function(x) read_tsv(file.path(filePathToData, x)))
  )

# store file names in a list
filePathToData <- '/Users/acadiak/Desktop/CGI/speciesPairs/CGIcentric/Noonan_summaryFiles/'
files <- dir(path = filePathToData, pattern = "*.txt")

# read in as nested dataframe
dataNoonan <- tibble(filename = files) %>%
  mutate(file_contents = map(filename,
                             function(x) read_tsv(file.path(filePathToData, x)))
  )

# combine Roller & Noonan tables
data <- full_join(dataRoller, dataNoonan)

# unnest to generate full table, but with file name as a column
unnested <- unnest(data, cols = c(file_contents))

# expand filename into full info
preFilter <- unnested %>% 
  # remove .txt from filename
  mutate(across('filename', str_replace, '.txt', '')) %>%
  # separate file names into variables
  separate(filename, into = c('speciesA', 'speciesB', 'tissue', 'mark', 'timePoint'))
preFilter$speciesPair <- paste(preFilter$speciesA, preFilter$speciesB, sep = '_')

table <- preFilter

# Rename species from genome abbreviations into names and rename age levels
ageLevels <- c('OlderThanAmniota', 'Amniota', 'Mammalia',
               'Theria', 'Eutheria', 'Primate', 'Ape', 'Human', 'None')

table <- table %>%
  mutate(speciesA = factor(speciesA, 
                           levels = c('hg19', 'rheMac10', 'calJac4', 'mm39', 'rn7', 'susScr11', 'canFam6', 'felCat9'), 
                           labels = c('Human', 'Rhesus', 'Marmoset', 'Mouse', 'Rat', 'Pig', 'Dog', 'Cat')),
         speciesB = factor(speciesB,
                           levels = c('rheMac2', 'calJac4', 'mm39', 'mm9', 'rn7', 'susScr11', 'canFam6', 'felCat9', 'equCab3'), 
                           labels = c('Rhesus', 'Marmoset', 'Mouse', 'Mouse', 'Rat', 'Pig', 'Dog', 'Cat', 'Horse')),
         OldestSegment_CGI = factor(OldestSegment_CGI,
                                    levels = rev(ageLevels), 
                                    labels = rev(c('Older Than Amniota', ageLevels[2:length(ageLevels)]))),
         OldestSegment_Peak = factor(OldestSegment_Peak,
                                     levels = rev(ageLevels), 
                                     labels = rev(c('Older Than Amniota', ageLevels[2:length(ageLevels)])))
  )

# rename speciesPair column now that they're renamed as species instead of genomes
table$speciesPair <- paste(table$speciesA, table$speciesB, sep = '\n')

# replace ages younger than Eutheria with "None", if >= 1 species in the pair is not a primate or ape/human
nonPrimates <- c('Mouse', 'Rat', 'Pig', 'Dog', 'Cat', 'Horse')
nonApe_Human <- c('Rhesus', 'Marmoset', 'Mouse', 'Rat', 'Pig', 'Dog', 'Cat', 'Horse')

# for CGI
table$OldestSegment_CGI_backup <- table$OldestSegment_CGI
table[(table$speciesA %in% nonPrimates | table$speciesB %in% nonPrimates) & table$OldestSegment_CGI=='Primate',]$OldestSegment_CGI <- "None"
table[(table$speciesA %in% nonPrimates | table$speciesB %in% nonPrimates) & table$OldestSegment_CGI=='Ape',]$OldestSegment_CGI <- "None"
table[(table$speciesA %in% nonPrimates | table$speciesB %in% nonPrimates) & table$OldestSegment_CGI=='Human',]$OldestSegment_CGI <- "None"

table[(table$speciesA %in% nonApe_Human | table$speciesB %in% nonApe_Human) & table$OldestSegment_CGI=='Ape',]$OldestSegment_CGI <- "None"
table[(table$speciesA %in% nonApe_Human | table$speciesB %in% nonApe_Human) & table$OldestSegment_CGI=='Human',]$OldestSegment_CGI <- "None"

# order factors
table$OldestSegment_CGI <- factor(table$OldestSegment_CGI,
                                  levels = rev(c('Older Than Amniota', 'Amniota', 'Mammalia',
                                                 'Theria', 'Eutheria', 'Primate', 'None')),
                                  labels = rev(c('Older Than Amniota', 'Amniota', 'Mammalia',
                                                 'Theria', 'Eutheria', 'Primate', 'Unknown')))

# add additional columns
table$CGI_summary <- '0'
table[table$CGI_A>0&table$CGI_B>0,]$CGI_summary <- 'Shared'
table[table$CGI_A>0&table$CGI_B==0,]$CGI_summary <- 'A-only'
table[table$CGI_A==0&table$CGI_B>0,]$CGI_summary <- 'B-only'
table$Peak_summary <- '0'
table[table$Peak_A>0&table$Peak_B>0,]$Peak_summary <- 'Shared'
table[table$Peak_A>0&table$Peak_B==0,]$Peak_summary <- 'A-only'
table[table$Peak_A==0&table$Peak_B>0,]$Peak_summary <- 'B-only'

# factor
table$CGI_summary <- factor(table$CGI_summary,
                            levels = c('A-only','B-only','Shared'))
table$Peak_summary <- factor(table$Peak_summary,
                            levels = c('A-only','B-only','Shared','0'))
table$mark <- factor(table$mark,
                     levels = c('H3K4me3', 'H3K27ac', 'H3K4me1', 'ac', 'me2'))
table$tissue <- factor(table$tissue,
                       levels = c('brain', 'liver', 'muscle', 'testis'),
                       labels = c('Brain', 'Liver', 'Muscle', 'Testis'))

# add column describing both CGI and Peak species-specificity
table$category <- paste(table$CGI_summary, table$Peak_summary, sep = '_')

# define function for getting 90% confidence interval, 25% quantile, median, 75% quantile
quantiles_90 <- function(x) {
  r <- quantile(x, probs=c(0.10, 0.25, 0.5, 0.75, 0.90))
  names(r) <- c("ymin", "lower", "middle", "upper", "ymax")
  r
}

# write this table to a text file since it takes so long to load it from the data
write_delim(table, 
            file = '/Users/acadiak/Desktop/CGI/R_scripts/Mar_2023/tables/table_peakStatComparisons.txt', 
            col_names = T,
            delim = '\t')

# read table back in
table <- read_tsv('/Users/acadiak/Desktop/CGI/R_scripts/Mar_2023/tables/table_peakStatComparisons.txt',
                  col_names = T) %>% tibble()

# factor - use this only if you read the table back in from the saved version
table$CGI_summary <- factor(table$CGI_summary,
                            levels = c('A-only','B-only','Shared'))
table$Peak_summary <- factor(table$Peak_summary,
                             levels = c('A-only','B-only','Shared','0'))
table$mark <- factor(table$mark,
                     levels = c('H3K4me3', 'H3K27ac', 'H3K4me1', 'ac', 'me2'))
table$tissue <- factor(table$tissue,
                       levels = c('Brain', 'Liver', 'Muscle', 'Testis'),
                       labels = c('Brain', 'Liver', 'Muscle', 'Testis'))

# store rhesus-mouse as a subset of the full table for memory
rhesusMouse <- table %>%
  filter(speciesPair == 'Rhesus\nMouse')

##### get info on phastCons for text
rhesusMouse$log2maxLOD <- log2(1 + rhesusMouse$MaxPhastConsScore_CGI)
shared <- rhesusMouse %>% filter(category == 'Shared_Shared', tissue == 'Brain', mark == 'H3K27ac')
A_only <- rhesusMouse %>% filter(category == 'A-only_A-only', tissue == 'Brain', mark == 'H3K27ac')
B_only <- rhesusMouse %>% filter(category == 'B-only_B-only', tissue == 'Brain', mark == 'H3K27ac')

round(quantile(shared$log2maxLOD, probs=c(0.10, 0.25, 0.5, 0.75, 0.90)), 2)
# 10%   25%   50%   75%   90% 
# 5.50  6.25  7.28  9.41 10.58 
round(quantile(A_only$log2maxLOD, probs=c(0.10, 0.25, 0.5, 0.75, 0.90)), 2)
# 10%  25%  50%  75%  90% 
# 0.00 3.81 5.32 7.39 9.09 
round(quantile(B_only$log2maxLOD, probs=c(0.10, 0.25, 0.5, 0.75, 0.90)), 2)
# 10%  25%  50%  75%  90% 
# 0.00 4.25 5.83 7.47 9.22 

# percentage that overlap ANY phastCons element
shared$anyPhastCons <- 'none'
A_only$anyPhastCons <- 'none'
B_only$anyPhastCons <- 'none'

shared[shared$log2maxLOD>0,]$anyPhastCons <- 'atLeastOne'
A_only[A_only$log2maxLOD>0,]$anyPhastCons <- 'atLeastOne'
B_only[B_only$log2maxLOD>0,]$anyPhastCons <- 'atLeastOne'

#                                 phastCons?    yes   no 
shared %>% group_by(anyPhastCons) %>% count() # 90    5      round(90/(90+5)*100, 2) = 94.74
A_only %>% group_by(anyPhastCons) %>% count() # 318   103    round(318/(318+103)*100, 2) = 75.53
B_only %>% group_by(anyPhastCons) %>% count() # 410   80     round(410/(410+80)*100, 2) = 83.67

#######################################################
# phastCons supp with many pairs phastCons comparison #
#######################################################

# Choose a subset of species pairs - 6
speciesPairList <- c('Rhesus\nMouse', 'Rat\nDog', 'Marmoset\nPig', 'Mouse\nRat', 'Pig\nCat', 'Dog\nHorse')

# wilcoxon test
wilcox.phastCons.maxLOD.A.vs.Shared <- table %>%
  filter(category %in% c('A-only_A-only', 'Shared_Shared'),
         speciesPair %in% speciesPairList) %>%
  group_by(speciesPair, tissue, mark) %>%
  do(w = wilcox.test(MaxPhastConsScore_CGI~category, data=., paired = F)) %>%
  summarize(speciesPair, tissue, mark, Wilcox = w$p.value)
wilcox.phastCons.maxLOD.A.vs.Shared$test <- 'A'

wilcox.phastCons.maxLOD.B.vs.Shared <- table %>%
  filter(category %in% c('B-only_B-only', 'Shared_Shared'),
         speciesPair %in% speciesPairList) %>%
  group_by(speciesPair, tissue, mark) %>%
  do(w = wilcox.test(MaxPhastConsScore_CGI~category, data=., paired = F)) %>%
  summarize(speciesPair, tissue, mark, Wilcox = w$p.value)
wilcox.phastCons.maxLOD.B.vs.Shared$test <- 'B'

# combine two tables and adjust p values
wilcox.phastCons.maxLOD <- full_join(wilcox.phastCons.maxLOD.A.vs.Shared, wilcox.phastCons.maxLOD.B.vs.Shared)
wilcox.phastCons.maxLOD$Wilcox.adj <- p.adjust(wilcox.phastCons.maxLOD$Wilcox, method = 'BH')

# add columns for plotting with geom_segment and geom_text
wilcox.phastCons.maxLOD$x <- 0
wilcox.phastCons.maxLOD[wilcox.phastCons.maxLOD$test=='A',]$x <- 1.2
wilcox.phastCons.maxLOD[wilcox.phastCons.maxLOD$test=='B',]$x <- 2.2
wilcox.phastCons.maxLOD$xend <- 2.8
wilcox.phastCons.maxLOD$y <- 0
wilcox.phastCons.maxLOD[wilcox.phastCons.maxLOD$test=='A',]$y <- 13.5
wilcox.phastCons.maxLOD[wilcox.phastCons.maxLOD$test=='B',]$y <- 11.5
wilcox.phastCons.maxLOD$yend <- wilcox.phastCons.maxLOD$y
wilcox.phastCons.maxLOD$stars <- ''
wilcox.phastCons.maxLOD[wilcox.phastCons.maxLOD$Wilcox.adj<0.05,]$stars <- '*'
wilcox.phastCons.maxLOD$stars.x <- (wilcox.phastCons.maxLOD$x + wilcox.phastCons.maxLOD$xend) / 2

# remove rows without star
wilcox.phastCons.maxLOD <- wilcox.phastCons.maxLOD[wilcox.phastCons.maxLOD$stars == '*',]

# store function for making plots
function_phastCons <- function(table, inMark, inColor) {
  x <- table %>%
  filter(speciesPair %in% speciesPairList,
         mark == inMark,
         category %in% c('A-only_A-only', 'B-only_B-only', 'Shared_Shared')) %>%
  ggplot(aes(x = category, y = log2(1 + MaxPhastConsScore_CGI), fill = inColor)) +
  facet_grid(tissue ~ fct_relevel(speciesPair, speciesPairList)) +
  stat_summary(fun.data = quantiles_90, 
               geom='boxplot', 
               position=position_dodge2(preserve = 'single'), width = 0.8, lwd = 0.25) +
  scale_fill_manual(values = inColor) +
  theme_bw() +
  theme(panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        strip.background = element_blank(),
        axis.text = element_text(size = 12),
        axis.text.x = element_text(angle = 0),
        axis.title = element_text(size = 12),
        legend.position = 'none',
        strip.text.x = element_text(size = 12)) +
  labs(y = 'log2(1 + Max phastCons\nLOD Score in oCGI',
       x = '') +
  coord_cartesian(ylim = c(0, 15)) +
  geom_segment(data = wilcox.phastCons.maxLOD %>% 
                 filter(speciesPair %in% speciesPairList,
                        mark == inMark),
               aes(x = x, xend = xend, 
                   y = y, yend = yend), 
               inherit.aes = F,
               lwd = 0.5) +
  geom_text(data = wilcox.phastCons.maxLOD %>% 
              filter(speciesPair %in% speciesPairList,
                     mark == inMark),
            aes(x = stars.x , y = y + 0.1, label = stars), inherit.aes = F, size = 4)
  return(x)
}

H3K4me3_phastCons <- function_phastCons(table, 'H3K4me3', 'tan1')
H3K27ac_phastCons <- function_phastCons(table, 'H3K27ac', 'darkseagreen3')
H3K4me1_phastCons <- function_phastCons(table, 'H3K4me1', 'mediumpurple2')

phastCons <- H3K4me3_phastCons + H3K27ac_phastCons + H3K4me1_phastCons + plot_layout(nrow = 3)
ggsave('FigS24_phastCons_manyPairs.pdf', phastCons, height = 3000, width = 1800, units = 'px')


# get examples for Fig 3 by restricting to Rhesus\nMouse
phastCons_rhesusMouse <- table %>%
  filter(speciesPair %in% c('Rhesus\nMouse'),
         tissue == 'Brain',
         category %in% c('A-only_A-only', 'B-only_B-only', 'Shared_Shared')) %>%
  ggplot(aes(x = category, y = log2(1 + MaxPhastConsScore_CGI), fill = mark)) +
  facet_grid(tissue ~ mark) +
  stat_summary(fun.data = quantiles_90, 
               geom='boxplot', 
               position=position_dodge2(preserve = 'single'), width = 0.8, lwd = 0.25) +
  scale_fill_manual(values = c('tan1', 'darkseagreen3', 'mediumpurple2')) +
  theme_bw() +
  theme(panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        strip.background = element_blank(),
        axis.text = element_text(size = 12),
        axis.text.x = element_text(angle = 0),
        axis.title = element_text(size = 12),
        legend.position = 'none',
        strip.text.x = element_text(size = 12)) +
  labs(y = 'log2(1 + Max phastCons\nLOD Score in oCGI',
       x = '') +
  coord_cartesian(ylim = c(0, 15)) +
  geom_segment(data = wilcox.phastCons.maxLOD %>% 
                 filter(speciesPair %in% 'Rhesus\nMouse',
                        tissue == 'Brain'),
               aes(x = x, xend = xend, 
                   y = y, yend = yend), 
               inherit.aes = F,
               lwd = 0.5) +
  geom_text(data = wilcox.phastCons.maxLOD %>% 
              filter(speciesPair %in% 'Rhesus\nMouse',
                     tissue == 'Brain'),
            aes(x = stars.x , y = y + 0.1, label = stars), inherit.aes = F, size = 6)
ggsave('Fig3E_phastCons_rhesusMouse.pdf', phastCons_rhesusMouse, height = 570, width = 1100, units = 'px')

# get median max LOD for rhesus_mouse comparison, for putting into the text

rhesusMouse %>%
  filter(speciesPair %in% c('Rhesus\nMouse'),
         tissue == 'Brain',
         mark == 'H3K27ac',
         category %in% c('A-only_A-only')) %>%
  dplyr::select(category, MaxPhastConsScore_CGI) %>%
  group_by(category) %>%
  summary() # median = 39.0

rhesusMouse %>%
  filter(speciesPair %in% c('Rhesus\nMouse'),
         tissue == 'Brain',
         mark == 'H3K27ac',
         category %in% c('B-only_B-only')) %>%
  dplyr::select(category, MaxPhastConsScore_CGI) %>%
  group_by(category) %>%
  summary() # median = 56.0

rhesusMouse %>%
  filter(speciesPair %in% c('Rhesus\nMouse'),
         tissue == 'Brain',
         mark == 'H3K27ac',
         category %in% c('Shared_Shared')) %>%
  dplyr::select(category, MaxPhastConsScore_CGI) %>%
  group_by(category) %>%
  summary() # median = 154.0


###########################################
# age supp with many pairs age comparison #
###########################################

# Choose a subset of species pairs - 6
speciesPairList <- c('Rhesus\nMouse', 'Rat\nDog', 'Marmoset\nPig', 'Mouse\nRat', 'Pig\nCat', 'Dog\nHorse')

# store function for making plots
function_age <- function(table, inMark) {
  x <- table %>%
    filter(speciesPair %in% speciesPairList,
           mark == inMark,
           category %in% c('A-only_A-only', 'B-only_B-only', 'Shared_Shared')) %>%
    group_by(speciesPair, tissue, category, OldestSegment_CGI) %>%
    count() %>%
    ggplot(aes(x = category, y = n, fill = OldestSegment_CGI)) +
    facet_grid(tissue ~ fct_relevel(speciesPair, speciesPairList)) +
    geom_bar(position = 'fill', stat = 'identity') + 
    scale_y_continuous(labels=scales::percent,
                       expand = c(0.01, 0),
                       breaks = c(0, 0.5, 1.0)) + 
    scale_fill_manual(values = c('gray90', 'gray70', '#c994c7', '#df65b0', '#dd1c77','#980043')) +
    theme_bw() + 
    theme(panel.grid.minor.x = element_blank(),
          panel.grid.major.x = element_blank(),
          panel.grid.minor.y = element_blank(),
          panel.spacing.x = unit(0.5, 'lines'),
          panel.spacing.y = unit(0.8, 'lines'),
          strip.background = element_blank(),
          axis.text = element_text(size = 10),
          #axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 12),
          #axis.text.x = element_blank(),
          axis.title = element_text(size = 20)) +
    labs(y = 'Percent of oCGIs with \nOldest Segment of Each Age',
         x = '')
  return(x)
}

H3K4me3_age <- function_age(table, 'H3K4me3')
H3K27ac_age <- function_age(table, 'H3K27ac')
H3K4me1_age <- function_age(table, 'H3K4me1')

age <- H3K4me3_age + H3K27ac_age + H3K4me1_age + plot_layout(nrow = 3)
ggsave('FigS25_age_manyPairs.pdf', age, height = 3000, width = 2000, units = 'px')


# get stats for text - number Amniota & older
rhesusMouse$age_binary <- 0
rhesusMouse[rhesusMouse$OldestSegment_CGI=='Amniota'|rhesusMouse$OldestSegment_CGI=='Older Than Amniota',]$age_binary <- 1
rhesusMouse_amniotaOlder <- rhesusMouse %>%
  filter(speciesPair %in% c('Rhesus\nMouse'),
         tissue == 'Brain',
         mark == 'H3K27ac',
         category %in% c('A-only_A-only', 'B-only_B-only', 'Shared_Shared')) %>%
  group_by(category, age_binary) %>%
  dplyr::count() %>%
  pivot_wider(names_from = age_binary,
              values_from = n)
rhesusMouse_amniotaOlder$pct <- rhesusMouse_amniotaOlder$'1' / (rhesusMouse_amniotaOlder$'1' + rhesusMouse_amniotaOlder$'0') * 100

# category        `0`   `1`   pct
# A-only_A-only   333    88  20.9
# B-only_B-only   401    89  18.2
# Shared_Shared    67    28  29.5

