# 2/13/23

require(ggplot2)
require(tidyverse)
require(patchwork)
require(cowplot)
theme_set(theme_cowplot())

setwd('/Users/acadiak/Desktop/CGI/Figures/Mar_2023')

# import data from singleSpecies CGIcentric table
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

# reorder species
table_CGI$species <- fct_relevel(table_CGI$filename, rev(c('Human', 'Rhesus', 'Marmoset', 'Mouse','Rat', 'Pig', 'Dog', 'Cat', 'Horse')))

# Fig S40 - percent of CGIs with phastBias in each species
table_CGI$phastBias <- '0'
table_CGI[table_CGI$CGI_phastBiasNumber>0,]$phastBias <- '1'
table_CGI$phastBias <- fct_relevel(table_CGI$phastBias, c('1', '0'))

FigS40 <- table_CGI %>%
  group_by(species, phastBias) %>%
  dplyr::count() %>%
  ggplot(aes(species, n/1000, fill = phastBias)) +
  geom_col(width = 0.6) +
  coord_flip() + 
  theme_minimal() +
  labs(y = 'Intronic and Intergenic\nCGIs (in Thousands)',
       x = '') +
  theme(axis.line = element_line(color = 'black', size = 1),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_line(size = 1),
        text = element_text(size = 12),
        axis.text = element_text(size = 12),
        legend.position = 'none') +
  scale_fill_manual(values = c('gray30', 'gray70')) +
  scale_y_continuous(limits = c(0, 85), expand = c(0, 0))

ggsave('FigS40_phastBiasOverview.pdf', FigS40, height=1080, width=850, units='px')

# get high and low percentages for text
percentageTable <- table_CGI %>%
  group_by(species, phastBias) %>%
  dplyr::count() %>%
  pivot_wider(names_from = phastBias, values_from = n)

percentageTable$pct <- percentageTable$'1' / (percentageTable$'1' + percentageTable$'0') * 100

View(percentageTable %>% dplyr::select(species, pct))

# Human 2.846408
# Rhesus 11.108568
# Mouse 13.969459
# Marmoset 26.203531
# Rat 33.225113
# Horse 51.941419
# Pig 63.409534
# Dog 66.175966
# Cat 67.265120
