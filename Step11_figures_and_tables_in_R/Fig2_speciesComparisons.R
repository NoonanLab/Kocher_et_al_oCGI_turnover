# 1/11/22
# Purpose:
# 1) generate plots in Fig 2C-E comparing CpG number, phastCons, and age in A-only, B-only, and Shared oCGIs
# 2) generate additional plots for supplement

# uses length and CpG num from original species coordinates (not merged oCGI)
# see Step5 for generation of these files in SummaryFiles_noPromFilter/

require(ggplot2)
require(cowplot)
require(tidyverse)
require(RColorBrewer)
theme_set(theme_cowplot())

# work here
setwd('/Users/acadiak/Desktop/CGI/Figures/Mar_2023')

# store files in a large table - each species pair has its own file
filePathToData <- '/Users/acadiak/Desktop/CGI/speciesPairs/CGIcentric/SummaryFiles_noPromFilter/'

# store file names in a list
files <- dir(path = filePathToData, pattern = "*.txt")

# read in as nested dataframe
data <- tibble(filename = files) %>%
  mutate(file_contents = map(filename,
                             function(x) read_tsv(file.path(filePathToData, x)))
  )

# unnest to generate full table, but with file name as a column
unnested <- unnest(data, cols = c(file_contents))

# modify variable names using various tidy commands
preFilter <- unnested %>% 
  # remove .txt from filename
  mutate(across('filename', str_replace, '.txt', '')) %>%
  # separate file names into variables
  separate(filename, into = c('speciesA', 'speciesB', 'tissue', 'mark', 'timePoint'))

# save as table - no filtering steps are required
table <- preFilter

# add additional columns describing CGI species-specificity
table$speciesPair <- paste(table$speciesA, table$speciesB, sep = '_')
table$CGI_summary <- '0'
table[table$CGI_A>0&table$CGI_B>0,]$CGI_summary <- 'Shared'
table[table$CGI_A>0&table$CGI_B==0,]$CGI_summary <- 'A-only'
table[table$CGI_A==0&table$CGI_B>0,]$CGI_summary <- 'B-only'

# Remove unnecessary columns (all dealing with intersection with peaks) and remove redundant rheMac2_mm9 comparison
table <- table %>%
  select(-contains(c('Peak', 'RPM', 'RPKM'))) %>%
  filter(speciesPair != 'rheMac2_mm9')

# Rename species from genome abbreviations into names
table <- table %>%
  mutate(speciesA = factor(speciesA, 
                           levels = c('hg19', 'rheMac10', 'calJac4', 'mm39', 'rn7', 'susScr11', 'canFam6', 'felCat9'), 
                           labels = c('Human', 'Rhesus', 'Marmoset', 'Mouse', 'Rat', 'Pig', 'Dog', 'Cat')),
         speciesB = factor(speciesB,
                           levels = c('rheMac2', 'calJac4', 'mm39', 'mm9', 'rn7', 'susScr11', 'canFam6', 'felCat9', 'equCab3'), 
                           labels = c('Rhesus', 'Marmoset', 'Mouse', 'Mouse', 'Rat', 'Pig', 'Dog', 'Cat', 'Horse')))

# rename speciesPair column now that they're renamed as species instead of genomes
table$speciesPair <- paste(table$speciesA, table$speciesB, sep = '\n')

# replace ages younger than Eutheria with "None", if >= 1 species in the pair is not a primate or ape/human
nonPrimates <- c('Mouse', 'Rat', 'Pig', 'Dog', 'Cat', 'Horse')
table$OldestSegment_CGI_backup <- table$OldestSegment_CGI
table[(table$speciesA %in% nonPrimates | table$speciesB %in% nonPrimates) & table$OldestSegment_CGI=='Primate',]$OldestSegment_CGI <- "None"
table[(table$speciesA %in% nonPrimates | table$speciesB %in% nonPrimates) & table$OldestSegment_CGI=='Ape',]$OldestSegment_CGI <- "None"
table[(table$speciesA %in% nonPrimates | table$speciesB %in% nonPrimates) & table$OldestSegment_CGI=='Human',]$OldestSegment_CGI <- "None"

nonApe_Human <- c('Rhesus', 'Marmoset', 'Mouse', 'Rat', 'Pig', 'Dog', 'Cat', 'Horse')
table[(table$speciesA %in% nonApe_Human | table$speciesB %in% nonApe_Human) & table$OldestSegment_CGI=='Ape',]$OldestSegment_CGI <- "None"
table[(table$speciesA %in% nonApe_Human | table$speciesB %in% nonApe_Human) & table$OldestSegment_CGI=='Human',]$OldestSegment_CGI <- "None"

# order factors
table$OldestSegment_CGI <- factor(table$OldestSegment_CGI,
                                  levels = rev(c('OlderThanAmniota', 'Amniota', 'Mammalia',
                                                 'Theria', 'Eutheria', 'Primate', 'None')),
                                  labels = rev(c('Older Than Amniota', 'Amniota', 'Mammalia',
                                                 'Theria', 'Eutheria', 'Primate', 'Unknown')))
table$CGI_summary <- factor(table$CGI_summary,
                           levels = c('A-only','B-only','Shared'))


#### Generate main figure (Fig 2C-E) and supplementary figure (Fig S2.3) using only data from Rhesus-Mouse species pair

# get median and upper quartile for max LOD score for text
table$log2MaxLOD <- log2(1 + table$MaxPhastConsScore_CGI)
shared <- table %>%
  filter(speciesPair == 'Rhesus\nMouse', CGI_summary == 'Shared')
A_only <- table %>%
  filter(speciesPair == 'Rhesus\nMouse', CGI_summary == 'A-only')
B_only <- table %>%
  filter(speciesPair == 'Rhesus\nMouse', CGI_summary == 'B-only')

round(quantile(shared$log2MaxLOD, probs=c(0.10, 0.25, 0.5, 0.75, 0.90)), 2)
# 10%   25%   50%   75%   90% 
# 5.09  6.63  8.42 10.14 11.27 
round(quantile(A_only$log2MaxLOD, probs=c(0.10, 0.25, 0.5, 0.75, 0.90)), 2)
# 10%  25%  50%  75%  90% 
# 0.00 0.00 4.95 7.11 9.01 
round(quantile(B_only$log2MaxLOD, probs=c(0.10, 0.25, 0.5, 0.75, 0.90)), 2)
# 10%  25%  50%  75%  90% 
# 0.00 3.81 5.49 7.61 9.46 

# get percentage of oCGIs in each category with log2maxLOD > 0 (i.e., they overlap any phastCons element)
shared$anyPhastCons <- 'none'
A_only$anyPhastCons <- 'none'
B_only$anyPhastCons <- 'none'
shared[shared$log2MaxLOD>0,]$anyPhastCons <- 'atLeastOne'
A_only[A_only$log2MaxLOD>0,]$anyPhastCons <- 'atLeastOne'
B_only[B_only$log2MaxLOD>0,]$anyPhastCons <- 'atLeastOne'

# results                         phastCons?    yes     no      pct_yes
shared %>% group_by(anyPhastCons) %>% count() # 777 vs 39       round(777/(777+39)*100, 2) = 95.22
A_only %>% group_by(anyPhastCons) %>% count() # 2914 vs 1435    round(2914/(2914+1435)*100, 2) = 67.00
B_only %>% group_by(anyPhastCons) %>% count() # 3916 vs 1129    round(3916/(3916+1129)*100, 2) = 77.62

# count # with oldest segment amniota or older for text
table$ageBinary <- 0
table[table$OldestSegment_CGI=='Amniota'|table$OldestSegment_CGI=='Older Than Amniota',]$ageBinary <- 1

ageAmniotaOlder <- table %>%
  filter(speciesPair == 'Rhesus\nMouse') %>%
  group_by(CGI_summary, ageBinary) %>%
  dplyr::count() %>%
  pivot_wider(names_from = ageBinary,
              values_from = n)
ageAmniotaOlder$pct <- ageAmniotaOlder$'1' / (ageAmniotaOlder$'1' + ageAmniotaOlder$'0') * 100

# CGI_summary   `0`   `1`   pct
# A-only       3593   756  17.4
# B-only       4104   941  18.7
# Shared        474   342  41.9

# define function for getting 90% confidence interval, 25% quantile, median, 75% quantile
quantiles_90 <- function(x) {
  r <- quantile(x, probs=c(0.10, 0.25, 0.5, 0.75, 0.90))
  names(r) <- c("ymin", "lower", "middle", "upper", "ymax")
  r
}

## factor species pairs in order to match bar plots in Fig S2.2
speciesPairSubsetA <- c('rheMac10_mm39', 'calJac4_rn7', 'rheMac10_rn7', 'calJac4_mm39', 'rheMac10_calJac4', 'mm39_rn7',
                        'susScr11_felCat9', 'canFam6_equCab3', 'susScr11_canFam6', 'canFam6_felCat9', 'felCat9_equCab3', 'susScr11_equCab3')
speciesPairSubsetB <- c('rheMac10_susScr11', 'calJac4_canFam6', 'mm39_felCat9', 'rn7_equCab3', 'rheMac10_canFam6', 'calJac4_felCat9', 'mm39_equCab3', 'rn7_susScr11',
                        'rheMac10_felCat9', 'calJac4_equCab3', 'mm39_susScr11', 'rn7_canFam6', 'rheMac10_equCab3', 'calJac4_susScr11', 'mm39_canFam6', 'rn7_felCat9')
speciesPairSubsetC <- c('hg19_rheMac2','hg19_mm9')

# substitute names
speciesPairList <- c(speciesPairSubsetA, speciesPairSubsetB, speciesPairSubsetC)
speciesPairList <- sub('hg19', 'Human', speciesPairList)
speciesPairList <- sub('rheMac10', 'Rhesus', speciesPairList)
speciesPairList <- sub('rheMac2', 'Rhesus', speciesPairList)
speciesPairList <- sub('calJac4', 'Marmoset', speciesPairList)
speciesPairList <- sub('mm39', 'Mouse', speciesPairList)
speciesPairList <- sub('mm9', 'Mouse', speciesPairList)
speciesPairList <- sub('rn7', 'Rat', speciesPairList)
speciesPairList <- sub('susScr11', 'Pig', speciesPairList)
speciesPairList <- sub('canFam6', 'Dog', speciesPairList)
speciesPairList <- sub('felCat9', 'Cat', speciesPairList)
speciesPairList <- sub('equCab3', 'Horse', speciesPairList)
speciesPairList <- sub('_', '\n', speciesPairList)

# factor
table$speciesPair <- factor(table$speciesPair, levels = speciesPairList)

#############################################
### GENERATE FIG 2 AND RELATED SUPPLEMENT ###
#############################################

##### Generate plots with CpG number for full set of species pairs 
##### for supplementary figure (Fig S12)

# do wilcoxon test
wilcox.CpGnum.A <- table %>%
  filter(CGI_summary %in% c('A-only', 'Shared')) %>%
  group_by(speciesPair) %>%
  do(w = wilcox.test(NumCpGs_CGI_A~CGI_summary, data=., paired = F)) %>%
  summarize(speciesPair, Wilcox = w$p.value)
wilcox.CpGnum.A$test <- 'A'

wilcox.CpGnum.B <- table %>%
  filter(CGI_summary %in% c('B-only', 'Shared')) %>%
  group_by(speciesPair) %>%
  do(w = wilcox.test(NumCpGs_CGI_B~CGI_summary, data=., paired = F)) %>%
  summarize(speciesPair, Wilcox = w$p.value)
wilcox.CpGnum.B$test <- 'B'

# combine two tables and adjust p values
wilcox.CpGnum <- full_join(wilcox.CpGnum.A, wilcox.CpGnum.B)
wilcox.CpGnum$Wilcox.adj <- p.adjust(wilcox.CpGnum$Wilcox, method = 'BH')

# add columns for plotting with geom_segment and geom_text
wilcox.CpGnum$x <- 1
wilcox.CpGnum$xend <- 2
wilcox.CpGnum$y <- 45
wilcox.CpGnum$yend <- wilcox.CpGnum$y
wilcox.CpGnum$stars <- ''
wilcox.CpGnum[wilcox.CpGnum$Wilcox.adj<0.05,]$stars <- '*'
wilcox.CpGnum$stars.x <- (wilcox.CpGnum$x + wilcox.CpGnum$xend) / 2

# drop rows without a star to prevent line from being plotted
wilcox.CpGnum <- wilcox.CpGnum[wilcox.CpGnum$Wilcox.adj < 0.05,]

# CpG numbers
allSpecies_CpGs_A <- table %>%
  filter(CGI_summary %in% c('A-only', 'Shared')) %>%
  ggplot(aes(x = CGI_summary,
             y = NumCpGs_CGI_A,
             fill = CGI_summary)) + 
  facet_grid(. ~ speciesPair, scales = 'free') +
  stat_summary(fun.data = quantiles_90, 
               geom='boxplot', 
               position=position_dodge2(preserve = 'single'), lwd = 0.25) +
  scale_fill_manual(values = c('royalblue4', 'gray70')) +  
  theme_bw() + 
  theme(panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        strip.background = element_blank(),
        axis.text.y = element_text(size = 12),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 8),
        axis.title = element_text(size = 20),
        legend.position = 'none') +
  labs(y = 'CpGs in\nSpecies A',
       x = '') +
  geom_segment(data = wilcox.CpGnum[wilcox.CpGnum$test=='A',], 
               aes(x = x, xend = xend, 
                   y = y, yend = yend), 
               inherit.aes = F,
               lwd = 0.5) +
  geom_text(data = wilcox.CpGnum[wilcox.CpGnum$test=='A',],
            aes(x = stars.x , y = y + 0.1, label = stars), inherit.aes = F, size = 6)

allSpecies_CpGs_B <- table %>%
  filter(CGI_summary %in% c('B-only', 'Shared')) %>%
  ggplot(aes(x = CGI_summary,
             y = NumCpGs_CGI_B,
             fill = CGI_summary)) + 
  facet_grid(. ~ speciesPair, scales = 'free') +
  stat_summary(fun.data = quantiles_90, 
               geom='boxplot', 
               position=position_dodge2(preserve = 'single'), lwd = 0.25) +
  scale_fill_manual(values = c('lightskyblue', 'gray70')) +  
  theme_bw() + 
  theme(panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        strip.background = element_blank(),
        axis.text.y = element_text(size = 12),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 8),
        axis.title = element_text(size = 20),
        legend.position = 'none') +
  labs(y = 'CpGs in\nSpecies B',
       x = '') +
  geom_segment(data = wilcox.CpGnum[wilcox.CpGnum$test=='B',], 
               aes(x = x, xend = xend, 
                   y = y, yend = yend), 
               inherit.aes = F,
               lwd = 0.5) +
  geom_text(data = wilcox.CpGnum[wilcox.CpGnum$test=='B',],
            aes(x = stars.x , y = y + 0.1, label = stars), inherit.aes = F, size = 6)

allSpecies_CpGs = allSpecies_CpGs_A / allSpecies_CpGs_B

ggsave('FigS12_allSpecies_CpGs.pdf', allSpecies_CpGs, height = 1200, width = 6000, units = 'px', limitsize = F)


# get medians for text for rhesus vs mouse CpG numbers - these were not actually used
summary(table %>%
          filter(speciesPair == 'Rhesus\nMouse', CGI_summary == 'A-only') %>% 
          select(NumCpGs_CGI_A)) # 13.0
summary(table %>%
          filter(speciesPair == 'Rhesus\nMouse', CGI_summary == 'Shared') %>% 
          select(NumCpGs_CGI_A)) # 15.0
summary(table %>%
          filter(speciesPair == 'Rhesus\nMouse', CGI_summary == 'B-only') %>% 
          select(NumCpGs_CGI_B)) # 11.0
summary(table %>%
          filter(speciesPair == 'Rhesus\nMouse', CGI_summary == 'Shared') %>% 
          select(NumCpGs_CGI_B)) # 13.0


##### Generate plots with LENGTH for full set of species pairs 
##### for supplementary figure (Fig S13)

# do wilcoxon test
wilcox.length.A <- table %>%
  filter(CGI_summary %in% c('A-only', 'Shared')) %>%
  group_by(speciesPair) %>%
  do(w = wilcox.test(lengthOriginalCGI_A~CGI_summary, data=., paired = F)) %>%
  summarize(speciesPair, Wilcox = w$p.value)
wilcox.length.A$test <- 'A'

wilcox.length.B <- table %>%
  filter(CGI_summary %in% c('B-only', 'Shared')) %>%
  group_by(speciesPair) %>%
  do(w = wilcox.test(lengthOriginalCGI_B~CGI_summary, data=., paired = F)) %>%
  summarize(speciesPair, Wilcox = w$p.value)
wilcox.length.B$test <- 'B'

# combine two tables and adjust p values
wilcox.length <- full_join(wilcox.length.A, wilcox.length.B)
wilcox.length$Wilcox.adj <- p.adjust(wilcox.length$Wilcox, method = 'BH')

# add columns for plotting with geom_segment and geom_text
wilcox.length$x <- 1.1
wilcox.length$xend <- 1.9
wilcox.length$y <- 900
wilcox.length$yend <- wilcox.length$y
wilcox.length$stars <- ''
wilcox.length[wilcox.length$Wilcox.adj<0.05,]$stars <- '*'
wilcox.length$stars.x <- (wilcox.length$x + wilcox.length$xend) / 2

# drop rows without a star to prevent line from being plotted
wilcox.length <- wilcox.length[wilcox.length$Wilcox.adj < 0.05,]

# length
allSpecies_length_A <- table %>%
  filter(CGI_summary %in% c('A-only', 'Shared')) %>%
  ggplot(aes(x = CGI_summary,
             y = lengthOriginalCGI_A,
             fill = CGI_summary)) + 
  facet_grid(. ~ speciesPair, scales = 'free') +
  stat_summary(fun.data = quantiles_90, 
               geom='boxplot', 
               position=position_dodge2(preserve = 'single'), lwd = 0.25) +
  scale_fill_manual(values = c('royalblue4', 'gray70')) +  
  theme_bw() + 
  theme(panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        strip.background = element_blank(),
        axis.text.y = element_text(size = 12),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 8),
        axis.title = element_text(size = 20),
        legend.position = 'none') +
  labs(y = 'Length in\nSpecies A',
       x = '') +
  geom_segment(data = wilcox.length[wilcox.length$test=='A',], 
               aes(x = x, xend = xend, 
                   y = y, yend = yend), 
               inherit.aes = F,
               lwd = 0.5) +
  geom_text(data = wilcox.length[wilcox.length$test=='A',],
            aes(x = stars.x , y = y + 0.1, label = stars), inherit.aes = F, size = 6)

allSpecies_length_B <- table %>%
  filter(CGI_summary %in% c('B-only', 'Shared')) %>%
  ggplot(aes(x = CGI_summary,
             y = lengthOriginalCGI_B,
             fill = CGI_summary)) + 
  facet_grid(. ~ speciesPair, scales = 'free') +
  stat_summary(fun.data = quantiles_90, 
               geom='boxplot', 
               position=position_dodge2(preserve = 'single'), lwd = 0.25) +
  scale_fill_manual(values = c('lightskyblue', 'gray70')) +  
  theme_bw() + 
  theme(panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        strip.background = element_blank(),
        axis.text.y = element_text(size = 12),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 8),
        axis.title = element_text(size = 20),
        legend.position = 'none') +
  labs(y = 'Length in\nSpecies B',
       x = '') +
  geom_segment(data = wilcox.length[wilcox.length$test=='B',], 
               aes(x = x, xend = xend, 
                   y = y, yend = yend), 
               inherit.aes = F,
               lwd = 0.5) +
  geom_text(data = wilcox.length[wilcox.length$test=='B',],
            aes(x = stars.x , y = y + 0.1, label = stars), inherit.aes = F, size = 6)

allSpecies_length = allSpecies_length_A / allSpecies_length_B

ggsave('FigS13_allSpecies_length.pdf', allSpecies_length, height = 1200, width = 6000, units = 'px', limitsize = F)


##### Generate plots with phastCons max LOD score for full set of species pairs
##### for supplementary figure (Fig S14)

# do wilcoxon test
wilcox.phastCons.maxLOD.A.vs.Shared <- table %>%
  filter(CGI_summary %in% c('A-only', 'Shared')) %>%
  group_by(speciesPair) %>%
  do(w = wilcox.test(MaxPhastConsScore_CGI~CGI_summary, data=., paired = F)) %>%
  summarize(speciesPair, Wilcox = w$p.value)
wilcox.phastCons.maxLOD.A.vs.Shared$test <- 'A'
  
wilcox.phastCons.maxLOD.B.vs.Shared <- table %>%
  filter(CGI_summary %in% c('B-only', 'Shared')) %>%
  group_by(speciesPair) %>%
  do(w = wilcox.test(MaxPhastConsScore_CGI~CGI_summary, data=., paired = F)) %>%
  summarize(speciesPair, Wilcox = w$p.value)
wilcox.phastCons.maxLOD.B.vs.Shared$test <- 'B'

# combine two tables and adjust p values
wilcox.phastCons.maxLOD <- full_join(wilcox.phastCons.maxLOD.A.vs.Shared, wilcox.phastCons.maxLOD.B.vs.Shared)
wilcox.phastCons.maxLOD$Wilcox.adj <- p.adjust(wilcox.phastCons.maxLOD$Wilcox, method = 'BH')

# add columns for plotting with geom_segment and geom_text
wilcox.phastCons.maxLOD$x <- 0
wilcox.phastCons.maxLOD[wilcox.phastCons.maxLOD$test=='A',]$x <- 1
wilcox.phastCons.maxLOD[wilcox.phastCons.maxLOD$test=='B',]$x <- 2
wilcox.phastCons.maxLOD$xend <- 3
wilcox.phastCons.maxLOD$y <- 0
wilcox.phastCons.maxLOD[wilcox.phastCons.maxLOD$test=='A',]$y <- 14
wilcox.phastCons.maxLOD[wilcox.phastCons.maxLOD$test=='B',]$y <- 12
wilcox.phastCons.maxLOD$yend <- wilcox.phastCons.maxLOD$y
wilcox.phastCons.maxLOD$stars <- '*'
wilcox.phastCons.maxLOD$stars.x <- (wilcox.phastCons.maxLOD$x + wilcox.phastCons.maxLOD$xend) / 2
  
# plot all species pairs
allSpecies_phastCons_maxLOD <- table %>%
  ggplot(aes(x = CGI_summary,
             y = log2(1+MaxPhastConsScore_CGI),
             fill = CGI_summary)) + 
  facet_wrap(vars(speciesPair), ncol = 6) +
  stat_summary(fun.data = quantiles_90, 
               geom='boxplot', 
               position=position_dodge2(preserve = 'single'), lwd = 0.25) +
  scale_fill_manual(values = c('royalblue4', 'lightskyblue', 'gray70')) +  
  theme_bw() + 
  theme(panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        strip.background = element_blank(),
        axis.text = element_text(size = 12),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 12),
        axis.title = element_text(size = 12),
        legend.position = 'none',
        strip.text.x = element_text(size = 12)) +
  labs(y = 'log2(1 + Max phastCons\nLOD Score in oCGI',
       x = '') +
  coord_cartesian(ylim = c(0, 15)) +
  geom_segment(data = wilcox.phastCons.maxLOD, 
                            aes(x = x, xend = xend, 
                                y = y, yend = yend), 
                            inherit.aes = F,
                            lwd = 0.5) +
  geom_text(data = wilcox.phastCons.maxLOD,
            aes(x = stars.x , y = y + 0.1, label = stars), inherit.aes = F, size = 6)

ggsave('FigS14_allSpecies_maxLOD.pdf', allSpecies_phastCons_maxLOD, height = 2400, width = 1800, units = 'px')

# compare medians for rhesus vs mouse to report in the text
summary(table %>% filter(speciesPair == 'Rhesus\nMouse', CGI_summary == 'A-only') %>% select(MaxPhastConsScore_CGI)) # 30.0
summary(table %>% filter(speciesPair == 'Rhesus\nMouse', CGI_summary == 'B-only') %>% select(MaxPhastConsScore_CGI)) # 44.0
summary(table %>% filter(speciesPair == 'Rhesus\nMouse', CGI_summary == 'Shared') %>% select(MaxPhastConsScore_CGI)) # 341.0

##### Generate plots with CGI age for full set of species pairs
##### for supplementary figure (Fig S16)

# count categories for ADULT tissues (Roller)
ageCounts <- table %>%
  group_by(speciesPair, CGI_summary, OldestSegment_CGI) %>%
  count()

allSpecies_age <- ageCounts %>%
  ggplot(aes(x = CGI_summary, y = n, fill = OldestSegment_CGI)) +
  facet_wrap(vars(speciesPair), ncol = 6) +
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

ggsave('FigS16_allSpecies_age.pdf', allSpecies_age, height = 2400, width = 2100, units = 'px')



##########
# Fig 2C #
##########

# CpG NUMBERS for Fig 2C - using CpG numbers from CpG islands in original species, not from merged interval
rhesusMouse_CpGnumbers_rhesus <- table %>%
  filter(CGI_summary %in% c('A-only', 'Shared'),
         speciesA=='Rhesus',
         speciesB=='Mouse') %>%
  ggplot(aes(x = CGI_summary,
             y = NumCpGs_CGI_A,
             fill = CGI_summary)) + 
  stat_summary(fun.data = quantiles_90, 
               geom='boxplot', 
               position=position_dodge2(preserve = 'single'), lwd = 0.25) +
  scale_fill_manual(values = c('royalblue4', 'gray70')) +  
  theme_bw() + 
  theme(panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        strip.background = element_blank(),
        axis.text = element_text(size = 12),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
        axis.title = element_text(size = 12),
        legend.position = 'none') +
  coord_cartesian(ylim = c(8, 50)) +
  labs(y = 'CpGs in Rhesus',
       x = '') +
  geom_segment(data = wilcox.CpGnum[wilcox.CpGnum$test=='A',] %>% filter(speciesPair == 'Rhesus\nMouse'), 
               aes(x = x, xend = xend, 
                   y = y, yend = yend), 
               inherit.aes = F,
               lwd = 0.5) +
  geom_text(data = wilcox.CpGnum[wilcox.CpGnum$test=='A',] %>% filter(speciesPair == 'Rhesus\nMouse'),
            aes(x = stars.x , y = y + 0.1, label = stars), inherit.aes = F, size = 6)

rhesusMouse_CpGnumbers_mouse <- table %>%
  filter(CGI_summary %in% c('B-only', 'Shared'),
         speciesA=='Rhesus',
         speciesB=='Mouse') %>%
  ggplot(aes(x = CGI_summary,
             y = NumCpGs_CGI_B,
             fill = CGI_summary)) + 
  stat_summary(fun.data = quantiles_90, 
               geom='boxplot', 
               position=position_dodge2(preserve = 'single'), lwd = 0.25) +
  scale_fill_manual(values = c('lightskyblue', 'gray70')) +  
  theme_bw() + 
  theme(panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        strip.background = element_blank(),
        axis.text = element_text(size = 12),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
        axis.title = element_text(size = 12),
        legend.position = 'none') +
  coord_cartesian(ylim = c(8, 50)) +
  labs(y = 'CpGs in Mouse',
       x = '') +
  geom_segment(data = wilcox.CpGnum[wilcox.CpGnum$test=='B',] %>% filter(speciesPair == 'Rhesus\nMouse'), 
               aes(x = x, xend = xend, 
                   y = y, yend = yend), 
               inherit.aes = F,
               lwd = 0.5) +
  geom_text(data = wilcox.CpGnum[wilcox.CpGnum$test=='B',] %>% filter(speciesPair == 'Rhesus\nMouse'),
            aes(x = stars.x , y = y + 0.1, label = stars), inherit.aes = F, size = 6)

Fig2C_rhesusMouse_CpGs <- rhesusMouse_CpGnumbers_rhesus + rhesusMouse_CpGnumbers_mouse
ggsave('Fig2C_rhesusMouse_CpGs.pdf', Fig2C_rhesusMouse_CpGs, height = 800, width = 800, units = 'px')

##########
# Fig 2D #
##########

# compare phastCons max LOD (Fig 2D)
Fig2D_rhesusMouse.phastCons.maxLOD <- table %>%
  ggplot(aes(x = CGI_summary,
             y = log2(1+MaxPhastConsScore_CGI),
             fill = CGI_summary)) + 
  stat_summary(fun.data = quantiles_90, 
               geom='boxplot', 
               position=position_dodge2(preserve = 'single'), lwd = 0.25) +
  scale_fill_manual(values = c('royalblue4', 'lightskyblue', 'gray70')) +  
  theme_bw() + 
  theme(panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        strip.background = element_blank(),
        axis.text = element_text(size = 12),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 12),
        axis.title = element_text(size = 12),
        legend.position = 'none',
        strip.text.x = element_text(size = 12)) +
  labs(y = 'log2(1 + Max phastCons\nLOD Score in oCGI',
       x = '') +
  coord_cartesian(ylim = c(0, 15)) +
  geom_segment(data = wilcox.phastCons.maxLOD %>% filter(speciesPair == 'Rhesus\nMouse'), 
               aes(x = x, xend = xend, 
                   y = y, yend = yend), 
               inherit.aes = F,
               lwd = 0.5) +
  geom_text(data = wilcox.phastCons.maxLOD %>% filter(speciesPair == 'Rhesus\nMouse'),
            aes(x = stars.x , y = y + 0.1, label = stars), inherit.aes = F, size = 6)

ggsave('Fig2D_rhesusMouse_phastCons.maxLOD.pdf', Fig2D_rhesusMouse.phastCons.maxLOD, height = 750, width = 550, units = 'px')

##########
# Fig 2E #
##########

# compare age distributions
rhesusMouse.age <- table %>%
  filter(speciesA=='Rhesus',
         speciesB=='Mouse') %>%
  group_by(speciesPair, CGI_summary, OldestSegment_CGI) %>%
  count() %>%
  ggplot(aes(x = CGI_summary,
             y = n,
             fill = OldestSegment_CGI)) +
  geom_bar(position = 'fill', stat = 'identity') + 
  scale_y_continuous(labels=scales::percent,
                     expand = c(0.01, 0)) + 
  scale_fill_manual(values = c('gray90', 'gray70', '#c994c7', '#df65b0', '#dd1c77','#980043')) +
  theme_bw() + 
  theme(panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        strip.background = element_blank(),
        axis.text = element_text(size = 20),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
        axis.title = element_text(size = 20)) +
  labs(y = 'Percent of oCGIs with \nOldest Segment of Each Age',
       x = '')

ggsave('Fig2E_rhesusMouse_age.pdf', rhesusMouse.age, height = 1300, width = 1400, units = 'px')


###########
# Fig S15 #
###########

# DO WILCOXON TESTS for all supplement S15 plots
  
# SUM of phastCons LOD SCORES
wilcox.rhesusMouse.sumLOD.A <- table %>%
  filter(CGI_summary %in% c('A-only', 'Shared'),
         speciesPair == 'Rhesus\nMouse') %>%
  group_by(speciesPair) %>%
  do(w = wilcox.test(SumOfPhastConsScores_CGI~CGI_summary, data=., paired = F)) %>%
  summarize(speciesPair, Wilcox = w$p.value)
wilcox.rhesusMouse.sumLOD.A$test <- 'sumLOD.A'
  
wilcox.rhesusMouse.sumLOD.B <- table %>%
  filter(CGI_summary %in% c('B-only', 'Shared'),
         speciesPair == 'Rhesus\nMouse') %>%
  group_by(speciesPair) %>%
  do(w = wilcox.test(SumOfPhastConsScores_CGI~CGI_summary, data=., paired = F)) %>%
  summarize(speciesPair, Wilcox = w$p.value)
wilcox.rhesusMouse.sumLOD.B$test <- 'sumLOD.B'

# add column to table for percent of Bp in phastCons element
table$phastConsPercentBp <- table$PhastConsBases_CGI / table$lengthCGI_human

wilcox.rhesusMouse.percentBp.A <- table %>%
  filter(CGI_summary %in% c('A-only', 'Shared'),
         speciesPair == 'Rhesus\nMouse') %>%
  group_by(speciesPair) %>%
  do(w = wilcox.test(phastConsPercentBp~CGI_summary, data=., paired = F)) %>%
  summarize(speciesPair, Wilcox = w$p.value)
wilcox.rhesusMouse.percentBp.A$test <- 'percentBp.A'
  
wilcox.rhesusMouse.percentBp.B <- table %>%
  filter(CGI_summary %in% c('B-only', 'Shared'),
         speciesPair == 'Rhesus\nMouse') %>%
  group_by(speciesPair) %>%
  do(w = wilcox.test(phastConsPercentBp~CGI_summary, data=., paired = F)) %>%
  summarize(speciesPair, Wilcox = w$p.value)
wilcox.rhesusMouse.percentBp.B$test <- 'percentBp.B'
  
# combine all and adjust p values
#wilcox.rhesusMouseSupp <- full_join(wilcox.rhesusMouse.length.A, wilcox.rhesusMouse.length.B)
#wilcox.rhesusMouseSupp <- full_join(wilcox.rhesusMouseSupp, wilcox.rhesusMouse.sumLOD.A)
wilcox.rhesusMouseSupp <- full_join(wilcox.rhesusMouse.sumLOD.A, wilcox.rhesusMouse.sumLOD.B)
wilcox.rhesusMouseSupp <- full_join(wilcox.rhesusMouseSupp, wilcox.rhesusMouse.percentBp.A)
wilcox.rhesusMouseSupp <- full_join(wilcox.rhesusMouseSupp, wilcox.rhesusMouse.percentBp.B)
wilcox.rhesusMouseSupp$Wilcox.adj <- p.adjust(wilcox.rhesusMouseSupp$Wilcox, method = 'BH')
  
# add info for plotting
wilcox.rhesusMouseSupp$x <- 1.1
wilcox.rhesusMouseSupp[wilcox.rhesusMouseSupp$test=='sumLOD.B',]$x <- 2.1
wilcox.rhesusMouseSupp[wilcox.rhesusMouseSupp$test=='percentBp.B',]$x <- 2.1
wilcox.rhesusMouseSupp$xend <- 2.9
#wilcox.rhesusMouseSupp[wilcox.rhesusMouseSupp$test=='length.A'|wilcox.rhesusMouseSupp$test=='length.B',]$xend <- 2
wilcox.rhesusMouseSupp$y <- 0
#wilcox.rhesusMouseSupp[wilcox.rhesusMouseSupp$test=='length.A'|wilcox.rhesusMouseSupp$test=='length.B',]$y <- 680
wilcox.rhesusMouseSupp[wilcox.rhesusMouseSupp$test=='sumLOD.A',]$y <- 14
wilcox.rhesusMouseSupp[wilcox.rhesusMouseSupp$test=='sumLOD.B',]$y <- 12
wilcox.rhesusMouseSupp[wilcox.rhesusMouseSupp$test=='percentBp.A',]$y <- 0.57
wilcox.rhesusMouseSupp[wilcox.rhesusMouseSupp$test=='percentBp.B',]$y <- 0.48
wilcox.rhesusMouseSupp$yend <- wilcox.rhesusMouseSupp$y

wilcox.rhesusMouseSupp$stars <- '*'
wilcox.rhesusMouseSupp$stars.x <- (wilcox.rhesusMouseSupp$x + wilcox.rhesusMouseSupp$xend) / 2

# MAKE PLOTS

# SUM OF phastCons LOD

FigS15_rhesusMouse_sumLOD <- table %>%
  ggplot(aes(x = CGI_summary,
             y = log2(1+SumOfPhastConsScores_CGI),
             fill = CGI_summary)) + 
  stat_summary(fun.data = quantiles_90, 
               geom='boxplot', 
               position=position_dodge2(preserve = 'single'), lwd = 0.25) +
  scale_fill_manual(values = c('royalblue4', 'lightskyblue', 'gray70')) +  
  theme_bw() + 
  theme(panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        strip.background = element_blank(),
        axis.text = element_text(size = 12),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 12),
        axis.title = element_text(size = 12),
        legend.position = 'none',
        strip.text.x = element_text(size = 12)) +
  labs(y = 'log2(1 + Sum of phastCons\nLOD Scores in oCGI',
       x = '') +
  coord_cartesian(ylim = c(0, 15)) +
  geom_segment(data = wilcox.rhesusMouseSupp %>% filter(test == 'sumLOD.A' | test == 'sumLOD.B'), 
               aes(x = x, xend = xend, 
                   y = y, yend = yend), 
               inherit.aes = F,
               lwd = 0.5) +
  geom_text(data = wilcox.rhesusMouseSupp %>% filter(test == 'sumLOD.A' | test == 'sumLOD.B'),
            aes(x = stars.x , y = y + 0.1, label = stars), inherit.aes = F, size = 6)

ggsave('FigS15_rhesusMouse_phastCons.sumLOD.pdf', FigS15_rhesusMouse_sumLOD, height = 750, width = 550, units = 'px')


# percent of Bp with phastCons element
FigS15_rhesusMouse_percentBp <- table %>%
  ggplot(aes(x = CGI_summary,
             y = phastConsPercentBp,
             fill = CGI_summary)) + 
  stat_summary(fun.data = quantiles_90, 
               geom='boxplot', 
               position=position_dodge2(preserve = 'single'), lwd = 0.25) +
  scale_fill_manual(values = c('royalblue4', 'lightskyblue', 'gray70')) +  
  theme_bw() + 
  theme(panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        strip.background = element_blank(),
        axis.text = element_text(size = 12),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 12),
        axis.title = element_text(size = 12),
        legend.position = 'none',
        strip.text.x = element_text(size = 12)) +
  labs(y = 'Percent of oCGI Overlapping\na phastCons Element',
       x = '') +
  scale_y_continuous(labels = scales::percent) +
  geom_segment(data = wilcox.rhesusMouseSupp %>% filter(test == 'percentBp.A' | test == 'percentBp.B'), 
               aes(x = x, xend = xend, 
                   y = y, yend = yend), 
               inherit.aes = F,
               lwd = 0.5) +
  geom_text(data = wilcox.rhesusMouseSupp %>% filter(test == 'percentBp.A' | test == 'percentBp.B'),
            aes(x = stars.x , y = y + 0.005, label = stars), inherit.aes = F, size = 6)

ggsave('FigS15_rhesusMouse_phastCons.percentBp.pdf', FigS15_rhesusMouse_percentBp, height = 750, width = 600, units = 'px')

