# 1/29/23
# Purpose: make plots to show relationship between CGIs and HGEs
# Uses files made in HGE_analysis.sh

require(ggplot2)
require(tidyverse)
require(cowplot)
require(patchwork)
theme_set(theme_cowplot())

setwd('/Users/acadiak/Desktop/CGI/Figures/Mar_2023')

# read results of resampling test
data <- read_delim('/Users/acadiak/Desktop/CGI/EnhancerClasses/R_HGEs/HGEs_resamplingSummary.txt',
                   delim = '\t',
                   col_names = F)

# rename columns
colnames(data) <- c('tissue', 'mark', 'timePoint', 'CGI_pattern', 'observed', 'expected', 'p')

# adjust p values
data$padj <- p.adjust(data$p, method = 'BH')

# convert to percents
table.int <- data %>%
  group_by(tissue, mark, timePoint) %>%
  mutate(percent_obs = observed / sum(observed),
         percent_exp = expected / sum(expected))
table.int$FD <- table.int$percent_obs / table.int$percent_exp
table.int$tissue_mark <- paste(table.int$tissue, table.int$mark, sep = '_')

## HEATMAP
table.int$shape <- NA
table.int[table.int$padj<0.05&table.int$FD>1,]$shape <- 17
table.int[table.int$padj<0.05&table.int$FD<1,]$shape <- 6

# factor
table.int$CGI_pattern <- factor(table.int$CGI_pattern,
                            levels = c('H', 'R', 'M', 'HR', 'HM', 'RM', 'HRM', 'None'))
table.int$timePoint <- factor(table.int$timePoint,
                              levels = rev(c('CS16', 'CS23', 'F2F', 'F2O',
                                         'E33', 'E41', 'E44', 'E47')))

# add log2FD
table.int$log2FD <- log2(table.int$FD)

# plotting in function
function_HGE_heatmap <- function(table, inTissueMark, highColor, lowLimit, highLimit) {
  x <- table %>%
  filter(tissue_mark == inTissueMark) %>%
  mutate(log2FD = replace(log2FD, log2FD > highLimit, highLimit),
         log2FD = replace(log2FD, log2FD < lowLimit, lowLimit),
         log2FD = replace(log2FD, log2FD == '-Inf', lowLimit)) %>%
  ggplot(aes(x = CGI_pattern, y = timePoint, fill = log2FD, shape = shape)) +
  facet_grid(tissue_mark ~ .) +
  geom_tile(color = 'white', size = 1) +
  scale_fill_gradient2(high = highColor, low = 'gray30', mid = 'white',
                       limits = c(lowLimit, highLimit)) +
  theme_bw() + 
  theme(panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        strip.background = element_blank(),
        axis.text = element_text(size = 12),
        axis.text.x = element_text(size = 12),
        axis.title = element_text(size = 12),
        axis.line = element_blank(),
        #legend.position = 'none'
  ) +
  labs(y = '',
       x = '') +
  geom_point(aes(x = CGI_pattern, y = timePoint, shape = shape), size = 2) +
  scale_shape_identity()
  return(x)
}

brain_ac <- function_HGE_heatmap(table.int, 'brain_ac', 'palegreen4', -2.1, 1.1)
brain_me2 <- function_HGE_heatmap(table.int, 'brain_me2', 'darkgoldenrod1', -2.1, 1.1)
limb_ac <- function_HGE_heatmap(table.int, 'limb_ac', 'palegreen4', -2.1, 1.1)

# now Fig S35
HGE <- brain_ac + brain_me2 + limb_ac + plot_layout(nrow = 3)
ggsave('FigS29_fullHGE.pdf', HGE, height = 2200, width = 1700, units = 'px')

# output supplementary table
suppTable <- table.int %>%
  mutate(tissue = factor(tissue,
                         levels = c('brain', 'limb'),
                         labels = c('Cortex', 'Limb')),
         mark = factor(mark,
                       levels = c('ac', 'me2'),
                       labels = c('H3K27ac', 'H3K4me2')),
         timePoint = factor(timePoint,
                            levels = c('CS16', 'CS23', 'F2F', 'F2O', 'E33', 'E41', 'E44', 'E47'),
                            labels = c('7 p.c.w', '8.5 p.c.w.', '12 p.c.w. (Frontal)', '12 p.c.w. (Occipital)', 'E33', 'E41', 'E44', 'E47')),
         CGI_pattern = factor(CGI_pattern,
                              levels = c('H', 'R', 'M', 'HR', 'HM', 'RM', 'HRM', 'None'),
                              labels = c('Human-only', 'Rhesus-only', 'Mouse-only', 'Human and Rhesus', 
                                         'Human and Mouse', 'Rhesus and Mouse', 'Human, Rhesus, and Mouse', 'None'))) %>%
  select(-c(FD, tissue_mark, shape))
colnames(suppTable) <- c('Tissue', 'Histone Modification', 'Time Point', 'oCGI Species Pattern',
                         'Number of HGEs with oCGI Species Pattern', 'Mean Number of Resampled non-HGEs with oCGI Species Pattern',
                         'p', 'q', 'Percent of HGEs with oCGI Species Pattern', 'Mean Percent of Resampled non-HGEs with oCGI Species Pattern',
                         'Log2(Enrichment in HGEs compared to non-HGEs)')
suppTable$`Percent of HGEs with oCGI Species Pattern` <- suppTable$`Percent of HGEs with oCGI Species Pattern` * 100
suppTable$`Mean Percent of Resampled non-HGEs with oCGI Species Pattern` <- suppTable$`Mean Percent of Resampled non-HGEs with oCGI Species Pattern` * 100

# write Table S6
write_delim(suppTable %>% arrange('Tissue', 'Histone Modification', 'Time Point', 'oCGI Species Pattern'),
            'TableS6_HGE.txt',
            delim = '\t',
            quote = 'none')


# use this for making the main figure panel

# pivot to long table
table <- table.int %>%
  select(-c(observed, expected)) %>%
  pivot_longer(cols = c(percent_obs, percent_exp),
               names_to = 'category',
               values_to = 'percent')
# combine tissue + mark into tissue_mark
table$tissue_mark <- paste(table$tissue, table$mark, sep = '_')

# factor variables
table$CGI_pattern <- factor(table$CGI_pattern,
                            levels = c('H', 'RM', 'R', 'M', 'HR', 'HM', 'HRM', 'None'))
table$category <- factor(table$category,
                         levels = c('percent_obs', 'percent_exp'))

# plot
Fig4B_H <- table %>%
  filter(tissue_mark %in% c('brain_ac'),
         CGI_pattern %in% c('H'),
         timePoint == 'CS23') %>%
  ggplot(aes(x = category, y = percent, fill = category)) +
  facet_grid(CGI_pattern ~ tissue_mark, scales = 'free') +
  geom_bar(position = 'dodge', stat = 'identity', width = 0.8) +
  scale_fill_manual(values = c('gray30', 'gray80')) +  
  scale_y_continuous(labels = scales::percent,
                     limits = c(0, 0.12),
                     breaks = c(0, 0.06, 0.12)) +
  theme_bw() + 
  theme(panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        strip.background = element_blank(),
        axis.text = element_text(size = 12),
        axis.text.x = element_blank(),
        axis.title = element_text(size = 12),
        legend.position = 'none'
  ) +
  labs(y = '',
       x = '') +
  geom_segment(aes(x = 1, xend = 2,
                   y = 0.09, yend = 0.09)) +
  geom_text(aes(x = 1.5, y = 0.095, label = '*'), size = 6)

Fig4B_RM <- table %>%
  filter(tissue_mark %in% c('brain_ac'),
         CGI_pattern %in% c('RM'),
         timePoint == 'CS23') %>%
  ggplot(aes(x = category, y = percent, fill = category)) +
  facet_grid(CGI_pattern ~ tissue_mark, scales = 'free') +
  geom_bar(position = 'dodge', stat = 'identity', width = 0.8) +
  scale_fill_manual(values = c('gray30', 'gray80')) +  
  scale_y_continuous(labels = scales::percent,
                     limits = c(0, 0.02),
                     breaks = c(0, 0.01, 0.02)) +
  theme_bw() + 
  theme(panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        strip.background = element_blank(),
        axis.text = element_text(size = 12),
        axis.text.x = element_blank(),
        axis.title = element_text(size = 12),
        legend.position = 'none'
  ) +
  labs(y = '',
       x = '')

Fig4B_HR <- table %>%
  filter(tissue_mark %in% c('brain_ac'),
         CGI_pattern %in% c('HR'),
         timePoint == 'CS23') %>%
  ggplot(aes(x = category, y = percent, fill = category)) +
  facet_grid(CGI_pattern ~ tissue_mark, scales = 'free') +
  geom_bar(position = 'dodge', stat = 'identity', width = 0.8) +
  scale_fill_manual(values = c('gray30', 'gray80')) +  
  scale_y_continuous(labels = scales::percent,
                     limits = c(0, 0.20),
                     breaks = c(0, 0.10, 0.20)) +
  theme_bw() + 
  theme(panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        strip.background = element_blank(),
        axis.text = element_text(size = 12),
        axis.text.x = element_blank(),
        axis.title = element_text(size = 12),
        legend.position = 'none'
  ) +
  labs(y = '',
       x = '') +
  geom_segment(aes(x = 1, xend = 2,
                   y = 0.17, yend = 0.17)) +
  geom_text(aes(x = 1.5, y = 0.178, label = '*'), size = 6)

Fig4B <- Fig4B_H / Fig4B_RM / Fig4B_HR
ggsave('Fig4B_HGE_enrichment.pdf', Fig4B, height = 1400, width = 600, units = 'px')


##### generate distributions for the supplemental figure (Fig S34 (originally S28) panel D)
# download files from /home/ak2267/project/EnhancerClasses/HGE/permutation:
# brain_ac_CS23_observed_HGEs_CGI_patterns.txt
# brain_ac_CS23_resampled_HGEs_CGI_patterns.txt
# they are here: /Users/acadiak/Desktop/CGI/EnhancerClasses/R_HGEs

# read in to generate histograms
observedExample <- read.table('/Users/acadiak/Desktop/CGI/EnhancerClasses/R_HGEs/brain_ac_CS23_observed_HGEs_CGI_patterns.txt', header=T)
expectedExample <- read.table('/Users/acadiak/Desktop/CGI/EnhancerClasses/R_HGEs/brain_ac_CS23_resampled_HGEs_CGI_patterns.txt', header=T)

# add totals to calculate percentages
observedExample$total <- observedExample$H + observedExample$R + observedExample$M + observedExample$HR + 
  observedExample$HM + observedExample$RM + observedExample$HRM + observedExample$None
expectedExample$total <- expectedExample$H + expectedExample$R + expectedExample$M + expectedExample$HR + 
  expectedExample$HM + expectedExample$RM + expectedExample$HRM + expectedExample$None

# generate histogram
pdf(file = 'HGE_exampleHistogram.pdf', height = 3, width = 4)
hist(expectedExample$H / expectedExample$total * 100, 
     breaks = 20,
     xlim = c(2, 10),
     ylim = c(0, 3000),
     las = 1,
     xlab = '',
     main = '')
abline(v = mean(expectedExample$H / expectedExample$total * 100), col = 'black', lwd = 2)
abline(v = observedExample$H / observedExample$total * 100, col = 'darkseagreen4', lwd = 2)
dev.off()

mean(expectedExample$H / expectedExample$total * 100) # 5.297398
observedExample$H / observedExample$total * 100 # 7.314815



