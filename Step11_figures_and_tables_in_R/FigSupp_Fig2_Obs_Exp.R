# 3/7/2024
# Purpose: generate plots matching Fig S12 but for Obs/Exp CpG ratios
# instead of CpG number

# model this code after code in Fig2_speciesComparisons.R

require(ggplot2)
require(cowplot)
require(tidyverse)
require(patchwork)
theme_set(theme_cowplot())

# work here
setwd('/Users/acadiak/Desktop/Yale/!Writing/Kocher_CpG/GB_Revisions/Revised_code/Step12_revisions/Obs_Exp')

# store files in a large table - each species pair has its own file
# these files were generated during revisions by adding C and G counts
# so that I can calculate Obs/Exp ratios
# see ak20240307_differenceRequirement.sh
# the specific files I'm using have C and G counts from ORIGINAL CpGs! not reconciled
filePathToData <- '/Users/acadiak/Desktop/Yale/!Writing/Kocher_CpG/GB_Revisions/Revised_code/Step12_revisions/Difference_requirement/SummaryFiles_noPromFilter_OriginalGCcounts'

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

# make speciesPairList to get pairs in same order as previous supplementary plots
speciesPairSubsetA <- c('rheMac10_mm39', 'calJac4_rn7', 'rheMac10_rn7', 'calJac4_mm39', 'rheMac10_calJac4', 'mm39_rn7',
                        'susScr11_felCat9', 'canFam6_equCab3', 'susScr11_canFam6', 'canFam6_felCat9', 'felCat9_equCab3', 'susScr11_equCab3')
speciesPairSubsetB <- c('rheMac10_susScr11', 'calJac4_canFam6', 'mm39_felCat9', 'rn7_equCab3', 'rheMac10_canFam6', 'calJac4_felCat9', 'mm39_equCab3', 'rn7_susScr11',
                        'rheMac10_felCat9', 'calJac4_equCab3', 'mm39_susScr11', 'rn7_canFam6', 'rheMac10_equCab3', 'calJac4_susScr11', 'mm39_canFam6', 'rn7_felCat9')
speciesPairSubsetC <- c('hg19_rheMac2','hg19_mm9')
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

# refactor speciesPair column
table$speciesPair <- factor(table$speciesPair, levels = speciesPairList)

# remove columns having to do with age (OldestSegment) or phastCons or phastBias
table <- table %>%
  select(!c(starts_with('Oldest'), starts_with('phast'),
            PhastConsBases_CGI, MaxPhastConsScore_CGI, SumOfPhastConsScores_CGI,
            tissue, mark, timePoint))

# refactor CGI
table$CGI_summary <- factor(table$CGI_summary,
                            levels = c('A-only','B-only','Shared'))


##### Generate plots with CpG number for full set of species pairs 
##### for supplementary figure (Fig S12)

# calculate Obs/Exp for each CGI
table$OE_A <- (table$NumCpGs_CGI_A * table$lengthOriginalCGI_A) / (table$NumC_A * table$NumG_A)
table$OE_B <- (table$NumCpGs_CGI_B * table$lengthOriginalCGI_B) / (table$NumC_B * table$NumG_B)

# define function for getting 90% confidence interval, 25% quantile, median, 75% quantile
quantiles_90 <- function(x) {
  r <- quantile(x, probs=c(0.10, 0.25, 0.5, 0.75, 0.90))
  names(r) <- c("ymin", "lower", "middle", "upper", "ymax")
  r
}


# do wilcoxon test
wilcox.OE.A <- table %>%
  filter(CGI_summary %in% c('A-only', 'Shared')) %>%
  group_by(speciesPair) %>%
  do(w = wilcox.test(OE_A~CGI_summary, data=., paired = F)) %>%
  summarize(speciesPair, Wilcox = w$p.value)
wilcox.OE.A$test <- 'A'

wilcox.OE.B <- table %>%
  filter(CGI_summary %in% c('B-only', 'Shared')) %>%
  group_by(speciesPair) %>%
  do(w = wilcox.test(OE_B~CGI_summary, data=., paired = F)) %>%
  summarize(speciesPair, Wilcox = w$p.value)
wilcox.OE.B$test <- 'B'

# combine two tables and adjust p values
wilcox.OE <- full_join(wilcox.OE.A, wilcox.OE.B)
wilcox.OE$Wilcox.adj <- p.adjust(wilcox.OE$Wilcox, method = 'BH')

# add columns for plotting with geom_segment and geom_text
wilcox.OE$x <- 1.1
wilcox.OE$xend <- 1.9
wilcox.OE$y <- 0.76
wilcox.OE$yend <- wilcox.OE$y
wilcox.OE$stars <- ''
wilcox.OE[wilcox.OE$Wilcox.adj<0.05,]$stars <- '*'
wilcox.OE$stars.x <- (wilcox.OE$x + wilcox.OE$xend) / 2

# drop rows without a star to prevent line from being plotted
wilcox.OE <- wilcox.OE[wilcox.OE$Wilcox.adj < 0.05,]


# plot
allSpecies_OE_A <- table %>%
  filter(CGI_summary %in% c('A-only', 'Shared')) %>%
  ggplot(aes(x = CGI_summary,
             y = OE_A,
             fill = CGI_summary)) + 
  #ylim(0,3) +
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
  labs(y = 'Ratio of Observed / Expected CpGs\nin Species A',
       x = '') +
  geom_segment(data = wilcox.OE[wilcox.OE$test=='A',], 
               aes(x = x, xend = xend, 
                   y = y, yend = yend), 
               inherit.aes = F,
               lwd = 0.5) +
  geom_text(data = wilcox.OE[wilcox.OE$test=='A',],
            aes(x = stars.x , y = y + 0.01, label = stars), inherit.aes = F, size = 6)

allSpecies_OE_B <- table %>%
  filter(CGI_summary %in% c('B-only', 'Shared')) %>%
  ggplot(aes(x = CGI_summary,
             y = OE_B,
             fill = CGI_summary)) + 
  #ylim(0,3) +
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
  geom_segment(data = wilcox.OE[wilcox.OE$test=='B',], 
               aes(x = x, xend = xend, 
                   y = y, yend = yend), 
               inherit.aes = F,
               lwd = 0.5) +
  geom_text(data = wilcox.OE[wilcox.OE$test=='B',],
            aes(x = stars.x , y = y + 0.01, label = stars), inherit.aes = F, size = 6)

allSpecies_OE = allSpecies_OE_A / allSpecies_OE_B

ggsave('FigS17_allSpecies_OE.pdf', allSpecies_OE, height = 1200, width = 5500, units = 'px', limitsize = F)

