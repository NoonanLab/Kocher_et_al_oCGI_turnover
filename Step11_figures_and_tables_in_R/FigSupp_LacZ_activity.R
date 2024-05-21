# 4/21/23
# Purpose: analyze LacZ data from VISTA to see how predictive each mark is of activity (with and without CGIs)
# This version is updated to include comparisons to non-active VISTA elements

# MAKES FIG S3 (originally S2)

library(cowplot)
theme_set(theme_cowplot())
library(ggplot2)
library(tools) # to use toTitleCase

setwd('/Users/acadiak/Desktop/CGI/Figures/Apr_2023')

filePathToData <- '/Users/acadiak/Desktop/CGI/ENCODE_LacZ/summaryTables_VISTA-centric/'

# store file names in a list
files <- dir(path = filePathToData, pattern = "*.txt")

# read in as nested dataframe
data <- tibble(filename = files) %>%
  mutate(file_contents = map(filename,
                             function(x) read_tsv(file.path(filePathToData, x)))
  )

# unnest to generate full table, but with file name as a column
unnested <- unnest(data, cols = c(file_contents))

# separate filename
table <- unnested %>%
  separate(filename, into = c('timePoint', 'tissue', 'mark'),
           remove = T, sep = '_')

# store order for tissues and marks and factor in that order
tissues <- c('forebrain','midbrain','hindbrain','heart','limb')
marks <- c('H3K4me3','H3K27ac','H3K4me1','H3K4me2')

table$tissue <- factor(table$tissue, levels = tissues, labels = toTitleCase(tissues))
table$mark <- factor(table$mark, levels = marks)
table$CGI <- factor(table$CGI, levels = c(1, 0), labels = c('+', '-'))
table$peakBinary <- factor(table$peakBinary, levels = c(1, 0))

# paste together CGI and peakBinary (whether VISTA overlaps a peak)
table$label <- paste0(as.character(table$CGI), '_', as.character(table$peakBinary))
table$label <- factor(table$label,
                      levels = c('-_1', '-_0', '+_1', '+_0'),
                      labels = c('noCGI_peak', 'noCGI_bkgd', 'yesCGI_peak', 'yesCGI_bkgd'))

# do Fisher exact tests for each tissue x mark x CGI status - comparing peak same vs not same to non-peak same vs not same
# add new column converting 3-option LacZ_result column to 2-option isItSame column
table$isItSame <- 0
table[table$LacZ_result=='sameTissue',]$isItSame <- 1

# make count table
fisher.input <- table %>%
  group_by(tissue, mark, CGI, isItSame, peakBinary) %>%
  count()

# table for fisher testing:
#         peak    noPeak
# same
# notSame

# table for storing results of fisher test
resultsTable <- data.frame(tissue=character(),
                           mark=character(),
                           CGI=character(),
                           p=numeric())

# perform fisher tests and store output in table
for (queryTissue in c('Forebrain', 'Midbrain', 'Hindbrain', 'Heart', 'Limb')) {
  for (queryMark in c('H3K4me3', 'H3K27ac', 'H3K4me1', 'H3K4me2')) {
    for (queryCGI in c('-', '+')) {
      test.array <- fisher.input %>%
        filter(tissue == queryTissue, mark == queryMark, CGI == queryCGI)
      peak_same <- test.array %>% filter(peakBinary == 1, isItSame == 1)
      peak_notSame <- test.array %>% filter(peakBinary == 1, isItSame == 0)
      noPeak_same <- test.array %>% filter(peakBinary == 0, isItSame == 1)
      noPeak_notSame <- test.array %>% filter(peakBinary == 0, isItSame == 0)
      p = fisher.test(matrix(c(peak_same$n,
                               peak_notSame$n,
                               noPeak_same$n,
                               noPeak_notSame$n),
                      nrow = 2))$p.value
      resultsTable[nrow(resultsTable)+1,] <- c(queryTissue, queryMark, queryCGI, p)
    }
  }
}

# correct p values to q values
resultsTable$q <- p.adjust(resultsTable$p, method = 'BH')

# add coordinates for plotting
resultsTable$x <- 1
resultsTable$xend <- 2
resultsTable$y <- 0.75
resultsTable$yend <- resultsTable$y
resultsTable$stars <- ''
resultsTable[resultsTable$q < 0.05, ]$stars <- '*'
resultsTable$stars.x <- (resultsTable$xend + resultsTable$x) / 2
resultsTable <- resultsTable[resultsTable$stars == '*',]

# summarize as percentages
percentTable <- table %>%
  group_by(tissue, mark, label, LacZ_result) %>%
  count() %>%
  pivot_wider(names_from = LacZ_result, values_from = n)
percentTable$total <- percentTable$negative + percentTable$otherTissue + percentTable$sameTissue
percentTable$pct_same <- percentTable$sameTissue / (percentTable$total) * 100
percentTable$pct_other <- percentTable$otherTissue / (percentTable$total) * 100

percentTable_long <- percentTable %>% pivot_longer(cols = c('pct_same', 'pct_other'), names_to = 'percent')

# factor
percentTable_long$mark <- factor(percentTable_long$mark,
                                 levels = c('H3K27ac', 'H3K4me3', 'H3K4me2', 'H3K4me1'))
percentTable_long$tissue <- factor(percentTable_long$tissue,
                                   levels = c('Forebrain', 'Midbrain', 'Hindbrain', 'Heart', 'Limb'))

# plot for no oCGIs
noCGIs <- percentTable_long %>%
  filter(label %in% c('noCGI_bkgd', 'noCGI_peak')) %>%
  ggplot(aes(x = factor(label), y = value, fill = percent)) +
  facet_grid(factor(mark) ~ factor(tissue)) +
  geom_col(aes(x = label, y = value/100, fill = percent), color = 'gray30', lwd = 0.2) +
  theme_bw() +
  theme(panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        strip.background = element_blank(),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 12),
        strip.text = element_text(size = 12)) +
  scale_y_continuous(labels = scales::percent,
                     breaks = c(0, 0.5, 1.0),
                     expand = c(0.01, 0.01),
                     limits = c(0, 1.0)) +
  scale_fill_manual(values = c('white', 'dodgerblue4')) +
  labs(y = 'Percent of Tested Enhancers with LacZ Activity',
       x = '') +
  geom_segment(data = resultsTable %>% 
                 filter(CGI == '-'),
               aes(x = x, xend = xend, 
                   y = y, yend = yend), 
               inherit.aes = F,
               lwd = 0.5, color = 'dodgerblue4') +
  geom_text(data = resultsTable %>% 
              filter(CGI == '-'),
            aes(x = stars.x , y = y + 0.03, label = stars), color = 'dodgerblue4', inherit.aes = F, size = 6)

# plot for yes oCGIs - note that there's a "no data" value for Heart H3K4me2 and that star needs to be removed
yesCGIs <- percentTable_long %>%
  filter(label %in% c('yesCGI_bkgd', 'yesCGI_peak')) %>%
  ggplot(aes(x = factor(label), y = value, fill = percent)) +
  facet_grid(factor(mark) ~ factor(tissue)) +
  geom_col(aes(x = label, y = value/100, fill = percent), color = 'gray30', lwd = 0.2) +
  theme_bw() +
  theme(panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        strip.background = element_blank(),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 12),
        strip.text = element_text(size = 12)) +
  scale_y_continuous(labels = scales::percent,
                     breaks = c(0, 0.5, 1.0),
                     expand = c(0.01, 0.01),
                     limits = c(0, 1.0)) +
  scale_fill_manual(values = c('white', 'dodgerblue4')) +
  labs(y = 'Percent of Tested Enhancers with LacZ Activity',
       x = '') +
  geom_segment(data = resultsTable %>% 
                 filter(CGI == '+'),
               aes(x = x, xend = xend, 
                   y = y, yend = yend), 
               inherit.aes = F,
               lwd = 0.5, color = 'dodgerblue4') +
  geom_text(data = resultsTable %>% 
              filter(CGI == '+'),
            aes(x = stars.x , y = y + 0.03, label = stars), color = 'dodgerblue4', inherit.aes = F, size = 6)

# now Fig S3
LacZ <- noCGIs + yesCGIs + plot_layout(ncol = 2)
ggsave('FigS2_LacZ.pdf', LacZ, height = 1800, width = 2800, units = 'px')

#### output counts as a table, Table S1
# table for fisher testing:
#         peak    noPeak
# same
# notSame

# columns: Tissue HistoneModification CGI NotActive_withPeak NotActive_withoutPeak Active_withPeak Active_withoutPeak p q
fisher.input_wide <- fisher.input %>%
  pivot_wider(names_from = c(isItSame, peakBinary),
              values_from = n)
fisher.input_wide <- left_join(fisher.input_wide, resultsTable %>% select(tissue, mark, CGI, p, q),
                               by = c('tissue', 'mark', 'CGI'))
colnames(fisher.input_wide) <- c('Tissue','Histone_Modification','CGI',
                                 'NotActive_withPeak', 'NotActive_withoutPeak',
                                 'Active_withPeak', 'Active_withoutPeak',
                                 'p', 'q')

# correct p and q to NA for heart H3K4me2 with CGI
fisher.input_wide[fisher.input_wide$Tissue=='Heart'&fisher.input_wide$Histone_Modification=='H3K4me2'&fisher.input_wide$CGI=='+',]$p <- NA
fisher.input_wide[fisher.input_wide$Tissue=='Heart'&fisher.input_wide$Histone_Modification=='H3K4me2'&fisher.input_wide$CGI=='+',]$q <- NA

# convert CGI indicators from +/- to 1/0
fisher.input_wide$CGI <- factor(fisher.input_wide$CGI,
                                levels = c('+', '-'),
                                labels = c('1', '0'))

# write table to output for use in supplement
write_delim(fisher.input_wide,
            file = 'TableS1_LacZ.txt',
            delim = '\t',
            col_names = T)
