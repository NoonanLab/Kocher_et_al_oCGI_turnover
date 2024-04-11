# 2/3/23
# Copy number qPCR to verify the mouse line

# THIS WAS FORMERLY FIGURE 5 - NOW REVISED FIGURE 4

require(tidyverse)
require(ggplot2)
require(cowplot)
theme_set(theme_cowplot())

setwd('/Users/acadiak/Desktop/CGI/Figures/Jan_2023')

# import table from sheet "Plate 7 & 8 combined results" in "qPCR_lineValidation_hs754"
qPCR <- read_tsv('~/Desktop/hs754-Misc/Validation/Copy\ number/copyNumber.txt')

# get control values by sample
control <- qPCR %>%
  filter(Primers == 'control')

# paste control back to table
norm_qPCR <- left_join(qPCR %>% filter(Primers != 'control'),
                       control,
                       by = c('Geno', 'Individual'),
                       suffix = c('.region', '.control'))

# add dCt, FC, and propgate SD
norm_qPCR$dCt <- norm_qPCR$Ct.region - norm_qPCR$Ct.control
norm_qPCR$FD <- 2 ^ ( - norm_qPCR$dCt)
norm_qPCR$SD.FD <- sqrt((norm_qPCR$SD.region^2) + (norm_qPCR$SD.control^2))

# add ddCt by normalizing to average dCt of all the WT samples
normTo1 <- norm_qPCR %>%
  filter(Geno == 'WT', Individual == '1') %>%
  select(Primers.region, dCt)

norm_qPCR <- left_join(norm_qPCR, normTo1,
                       by = 'Primers.region',
                       suffix = c('.preNorm', 'WT1'))

# calculate ddCt by normalizing to WT 1 and convert ddCt to FD.ddCt
norm_qPCR$ddCt <- norm_qPCR$dCt.preNorm - norm_qPCR$dCtWT1
norm_qPCR$FD.ddCt <- 2 ^ (- norm_qPCR$ddCt)

# scale SD.FD by the ratio of FD / FD.ddCt
norm_qPCR$SD.FD.ddCt <- norm_qPCR$SD.FD * (norm_qPCR$FD / norm_qPCR$FD.ddCt)

# factor variables
norm_qPCR$Geno <- factor(norm_qPCR$Geno,
                           levels = c('WT', 'HUM'))
norm_qPCR$Primers.region <- factor(norm_qPCR$Primers.region,
                            levels = c('hum', '5prime', '3prime', 'outside'),
                            labels = c('Edited Region', "5' Region", "3' Region", 'Adjacent Unedited Region'))

# plot
copyNumber <- norm_qPCR %>%
  ggplot(aes(x = Geno, y = FD.ddCt, fill = factor(Individual))) + 
  facet_grid(. ~ Primers.region) +
  geom_col(position = 'dodge', color = 'gray30') +
  geom_errorbar(aes(x = Geno, 
                    ymin = FD.ddCt - SD.FD.ddCt, 
                    ymax = FD.ddCt + SD.FD.ddCt),
                position = position_dodge(0.9),
                width = 0.4) +
  coord_cartesian(ylim = c(0, 1.5)) +
  scale_fill_manual(values = rep('gray80', 3)) +
  theme_bw() +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        strip.background = element_blank(),
        axis.text = element_text(size = 12),
        strip.text = element_text(size = 12),
        axis.title = element_text(size = 12),
        legend.position = 'none') +
  geom_hline(yintercept = 1, color = 'gray30', lty = 2) +
  labs(x = '', y = 'Fold Change vs Control Region\n(2 ^ -ddCt, Normalized to WT 1)')

ggsave('FigS30D_copyNumber_qPCR.pdf', copyNumber, height = 700, width = 2200, units = 'px')
