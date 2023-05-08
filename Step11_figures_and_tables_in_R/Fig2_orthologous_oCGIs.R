# 10/9/22
# Purpose:
# 1) generate bar plot with oCGI percentages for a subset of species pairs (Fig 2B)
# 2) generate bar plot with oCGI percentages for all species pairs (Fig S11)
# 3) generate Table S2 (numbers and percentages for all species pairs)

# updated 1/9/23 to use full tables to get counts
# see Step5 for generation of these files in SummaryFiles_noPromFilter/

# upset plots made with help from KY

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
preFilter$speciesPair <- paste(preFilter$speciesA, preFilter$speciesB, sep = '_')

# save as table - no filtering steps are required
table <- preFilter

# add additional columns
table$speciesPair <- paste(table$speciesA, table$speciesB, sep = '_')
table$CGI_summary <- '0'
table[table$CGI_A>0&table$CGI_B>0,]$CGI_summary <- 'both'
table[table$CGI_A>0&table$CGI_B==0,]$CGI_summary <- 'A'
table[table$CGI_A==0&table$CGI_B>0,]$CGI_summary <- 'B'

# remove unnecessary columns (all dealing with intersection with peaks)
table <- table %>%
  select(-contains(c('Peak', 'RPM', 'RPKM'))) %>%
  filter(speciesPair != 'rheMac2_mm9')


##########
# Fig 2B #
##########

##### Generate main figure bar plots (Fig 2B)

# store desired species pairs
speciesPairSubsetA <- c('rheMac10_mm39', 'calJac4_rn7', 'rheMac10_calJac4', 'mm39_rn7')
speciesPairSubsetB <- c('susScr11_canFam6', 'canFam6_felCat9', 'felCat9_equCab3', 'susScr11_felCat9', 'canFam6_equCab3')

# generate counts table
counts <- table %>%
  group_by(speciesPair, CGI_summary) %>%
  dplyr::count() %>%
  # make columns for speciesA and speciesB
  separate(speciesPair, into = c('speciesA', 'speciesB'), remove = F) %>%
  # convert counts to percents
  group_by(speciesPair) %>%
  mutate(pct = n/sum(n)) %>%
  ungroup

# add column defining subsets
counts$subset <- '0'
counts$subset[counts$speciesPair %in% speciesPairSubsetA] <- 'A'
counts$subset[counts$speciesPair %in% speciesPairSubsetB] <- 'B'

# remove NA subsets
counts <- counts[counts$subset!='0', ]

# set order of variables
counts$CGI_summary <- factor(counts$CGI_summary, levels = c('A','B','both'))
counts$speciesPair <- factor(counts$speciesPair, levels = c(speciesPairSubsetA, speciesPairSubsetB))

# make bar plot
barPlotMain <- counts %>%
  ggplot(aes(x = speciesPair, y = pct, fill = CGI_summary)) + 
  facet_grid(. ~ subset, scales = 'free', space = 'free') +
  geom_bar(stat = 'identity', color = 'black', width = 0.7, size = 1) +
  scale_fill_manual(values = c('white','gray20','gray50')) + 
  scale_y_continuous(expand = c(0, 0.0),
                     labels = scales::percent) +
  #theme_bw() +
  theme(strip.background = element_blank(), 
        #panel.grid.minor.y = element_blank(),
        strip.text = element_blank(),
        #axis.line.x = element_blank(), 
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 28)) + # element_text(angle = 90)) to add labels
  labs(y = 'Percent of Orthologous CGI Locations', fill = '', 
       x = '')

# define species subsets
speciesSubsetA <- c('rheMac10', 'calJac4', 'mm39', 'rn7')
speciesSubsetB <- c('susScr11', 'canFam6', 'felCat9', 'equCab3')

# species pairs for upset plot
upsetPairs <- within(counts, {
  speciesA.index <- speciesB.index <- NA
  speciesA.index[subset == 'A'] <- as.numeric(factor(speciesA[subset == 'A'], levels = speciesSubsetA))
  speciesB.index[subset == 'A'] <- as.numeric(factor(speciesB[subset == 'A'], levels = speciesSubsetA))
  speciesA.index[subset == 'B'] <- as.numeric(factor(speciesA[subset == 'B'], levels = speciesSubsetB))
  speciesB.index[subset == 'B'] <- as.numeric(factor(speciesB[subset == 'B'], levels = speciesSubsetB))
})

# empty circles for upset plot
upsetGrid <- expand.grid(speciesPair = c(speciesPairSubsetA, speciesPairSubsetB), 
                         species = 1:4) # or 1:max(length(speciesSubsetA), length(speciesSubsetB))
upsetGrid$subset[upsetGrid$speciesPair %in% speciesPairSubsetA] <- 'A'
upsetGrid$subset[upsetGrid$speciesPair %in% speciesPairSubsetB] <- 'B'

# make upsetPlot
upsetPlotMain <- ggplot(counts, aes(x = speciesPair, xend = speciesPair)) +
  facet_grid(. ~ subset, scales = 'free', space = 'free') +
  # plot gray bars
  geom_hline(yintercept = c(1, 3), size = 15, color = 'gray95') +
  # plot grey circles
  geom_point(data = upsetGrid, mapping = aes(y = species), 
             pch = 21, size = 8, stroke = 1, color = 'gray50', fill = 'gray90') +
  # upset plot
  geom_segment(data = upsetPairs, 
               mapping = aes(y = speciesA.index, yend = speciesB.index),
               lwd = 2) +
  geom_point(data = upsetPairs, mapping = aes(y = speciesA.index),
             pch = 21, size = 8, stroke = 2, fill = 'white') +
  geom_point(data = upsetPairs, mapping = aes(y = speciesB.index),
             pch = 21, size = 8, stroke = 2, fill = 'black') +
  # scales
  scale_y_reverse(expand = c(0, 0.5)) +
  # theme elements
  theme_void() + 
  theme(strip.text = element_blank(),
        panel.grid.major.x = element_line(size = 0.5, color = 'lightgray'))

## combine plots
Fig2B <- barPlotMain / upsetPlotMain + plot_layout(heights = c(2, 1))
ggsave('Fig2B_barplot.pdf', Fig2B, height=2500, width=3000, units='px')

# get numbers for text
counts %>%
  filter(speciesPair %in% c('rheMac10_mm39', 'rheMac10_calJac4', 'rn7_equCab3'))

# rheMac10_mm39 rheMac10 mm39     A            4349 0.426 
# rheMac10_mm39 rheMac10 mm39     B            5045 0.494 
# rheMac10_mm39 rheMac10 mm39     both          816 0.0799
# rheMac10_calJac4 rheMac10 calJac4  A            9643 0.435 
# rheMac10_calJac4 rheMac10 calJac4  B            8284 0.374 
# rheMac10_calJac4 rheMac10 calJac4  both         4241 0.191 
# rn7_equCab3      rn7      equCab3  A            9263 0.366 
# rn7_equCab3      rn7      equCab3  B           14177 0.560 
# rn7_equCab3      rn7      equCab3  both         1889 0.0746

###########
# Fig S11 #
###########

############################################################
#### Generate supplementary bar plots with all species pairs

# store desired species pairs
speciesPairSubsetA <- c('rheMac10_mm39', 'calJac4_rn7', 'rheMac10_rn7', 'calJac4_mm39', 'rheMac10_calJac4', 'mm39_rn7',
                        'susScr11_felCat9', 'canFam6_equCab3', 'susScr11_canFam6', 'canFam6_felCat9', 'felCat9_equCab3', 'susScr11_equCab3')
speciesPairSubsetB <- c('rheMac10_susScr11', 'calJac4_canFam6', 'mm39_felCat9', 'rn7_equCab3', 'rheMac10_canFam6', 'calJac4_felCat9', 'mm39_equCab3', 'rn7_susScr11',
                        'rheMac10_felCat9', 'calJac4_equCab3', 'mm39_susScr11', 'rn7_canFam6', 'rheMac10_equCab3', 'calJac4_susScr11', 'mm39_canFam6', 'rn7_felCat9')
speciesPairSubsetC <- c('hg19_rheMac2','hg19_mm9')

# generate counts table
counts <- table %>%
  group_by(speciesPair, CGI_summary) %>%
  count() %>%
  # make columns for speciesA and speciesB
  separate(speciesPair, into = c('speciesA', 'speciesB'), remove = F) %>%
  # convert counts to percents
  group_by(speciesPair) %>%
  mutate(pct = n/sum(n)) %>%
  ungroup

# add column defining subsets
counts$subset <- '0'
counts$subset[counts$speciesPair %in% speciesPairSubsetA] <- 'A'
counts$subset[counts$speciesPair %in% speciesPairSubsetB] <- 'B'
counts$subset[counts$speciesPair %in% speciesPairSubsetC] <- 'C'

# remove NA subsets
counts <- counts[counts$subset!='0', ]

# set order of variables
counts$CGI_summary <- factor(counts$CGI_summary, levels = c('A','B','both'))
counts$speciesPair <- factor(counts$speciesPair, levels = c(speciesPairSubsetA, speciesPairSubsetB, speciesPairSubsetC))

# make bar plot
barPlotSupp <- counts %>%
  ggplot(aes(x = speciesPair, y = pct, fill = CGI_summary)) + 
  facet_grid(. ~ subset, scales = 'free', space = 'free') +
  geom_bar(stat = 'identity', color = 'black', width = 0.7, size = 1) +
  scale_fill_manual(values = c('white','gray20','gray50')) + 
  scale_y_continuous(expand = c(0, 0.0),
                     labels = scales::percent) +
  #theme_bw() +
  theme(strip.background = element_blank(), 
        #panel.grid.minor.y = element_blank(),
        strip.text = element_blank(),
        #axis.line.x = element_blank(), 
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 36)) + # element_text(angle = 90)) to add labels
  labs(y = 'Percent of Orthologous CGI Locations', fill = '', 
       x = '')

# define species subsets
speciesSubsetA <- c('rheMac10', 'calJac4', 'mm39', 'rn7', 'susScr11', 'canFam6', 'felCat9', 'equCab3')
speciesSubsetB <- c('rheMac10', 'calJac4', 'mm39', 'rn7', 'susScr11', 'canFam6', 'felCat9', 'equCab3')
speciesSubsetC <- c('hg19', 'rheMac2', 'mm9')

# species pairs for upset plot
upsetPairs <- within(counts, {
  speciesA.index <- speciesB.index <- NA
  speciesA.index[subset == 'A'] <- as.numeric(factor(speciesA[subset == 'A'], levels = speciesSubsetA))
  speciesB.index[subset == 'A'] <- as.numeric(factor(speciesB[subset == 'A'], levels = speciesSubsetA))
  speciesA.index[subset == 'B'] <- as.numeric(factor(speciesA[subset == 'B'], levels = speciesSubsetB))
  speciesB.index[subset == 'B'] <- as.numeric(factor(speciesB[subset == 'B'], levels = speciesSubsetB))
  speciesA.index[subset == 'C'] <- as.numeric(factor(speciesA[subset == 'C'], levels = speciesSubsetC))
  speciesB.index[subset == 'C'] <- as.numeric(factor(speciesB[subset == 'C'], levels = speciesSubsetC))
})

# empty circles for upset plot
upsetGrid <- expand.grid(speciesPair = c(speciesPairSubsetA, speciesPairSubsetB, speciesPairSubsetC), 
                         species = 1:8) # or 1:max(length(speciesSubsetA), length(speciesSubsetB))
upsetGrid$subset[upsetGrid$speciesPair %in% speciesPairSubsetA] <- 'A'
upsetGrid$subset[upsetGrid$speciesPair %in% speciesPairSubsetB] <- 'B'
upsetGrid$subset[upsetGrid$speciesPair %in% speciesPairSubsetC] <- 'C'

# make upsetPlot
upsetPlotSupp <- ggplot(counts, aes(x = speciesPair, xend = speciesPair)) +
  facet_grid(. ~ subset, scales = 'free', space = 'free') +
  # plot gray bars
  geom_hline(yintercept = c(1, 3, 5, 7), size = 15, color = 'gray95') +
  # plot grey circles
  geom_point(data = upsetGrid, mapping = aes(y = species), 
             pch = 21, size = 8, stroke = 1, color = 'gray50', fill = 'gray90') +
  # upset plot
  geom_segment(data = upsetPairs, 
               mapping = aes(y = speciesA.index, yend = speciesB.index),
               lwd = 2) +
  geom_point(data = upsetPairs, mapping = aes(y = speciesA.index),
             pch = 21, size = 8, stroke = 2, fill = 'white') +
  geom_point(data = upsetPairs, mapping = aes(y = speciesB.index),
             pch = 21, size = 8, stroke = 2, fill = 'black') +
  # scales
  scale_y_reverse(expand = c(0, 0.5)) +
  # theme elements
  theme_void() + 
  theme(strip.text = element_blank(),
        panel.grid.major.x = element_line(size = 0.5, color = 'lightgray'))

## combine plots
Fig2_suppBarPlots <- barPlotSupp / upsetPlotSupp + plot_layout(heights = c(1, 1))
ggsave('FigS11_suppBarPlots.pdf', Fig2_suppBarPlots, height=4000, width=6000, units='px')


############
# Table S2 #
############

#### Generate supplementary table with full counts

counts_forTable <- counts %>%
  select(speciesA, speciesB, CGI_summary, n, pct) %>%
  mutate(pct = pct * 100) %>%
  mutate_if(is.numeric, round, digits = 2) %>%
  mutate(speciesA = factor(speciesA, 
                           levels = c('hg19', 'rheMac10', 'calJac4', 'mm39', 'rn7', 'susScr11', 'canFam6', 'felCat9', 'equCab3'), 
                           labels = c('Human', 'Rhesus', 'Marmoset', 'Mouse', 'Rat', 'Pig', 'Dog', 'Cat', 'Horse')),
         speciesB = factor(speciesB,
                           levels = c('rheMac2', 'calJac4', 'mm39', 'mm9', 'rn7', 'susScr11', 'canFam6', 'felCat9', 'equCab3'), 
                           labels = c('Rhesus', 'Marmoset', 'Mouse', 'Mouse', 'Rat', 'Pig', 'Dog', 'Cat', 'Horse'))) %>%
  pivot_wider(names_from = CGI_summary,
              values_from = c(n, pct)) %>%
  rename('Species_A' = 'speciesA',
         'Species_B' = 'speciesB',
         'A_only_oCGIs' = 'n_A',
         'B_only_oCGIs' = 'n_B',
         'Shared_oCGIs' = 'n_both',
         'A_only_percent' = 'pct_A',
         'B_only_percent' = 'pct_B',
         'Shared_percent' = 'pct_both') %>%
  arrange(factor(Species_A, levels = c('Human', 'Rhesus', 'Marmoset', 'Mouse', 'Rat', 'Pig', 'Dog', 'Cat', 'Horse')),
          factor(Species_B, levels = c('Human', 'Rhesus', 'Marmoset', 'Mouse', 'Rat', 'Pig', 'Dog', 'Cat', 'Horse')))

# write tab-delimited table (Table S2)
write_delim(counts_forTable,
            'TableS2_orthologous_oCGI_counts.txt',
            delim = '\t',
            quote = 'none')
