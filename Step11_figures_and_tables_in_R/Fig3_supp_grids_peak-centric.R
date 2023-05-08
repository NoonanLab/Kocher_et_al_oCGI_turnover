# 2/14/23
# PEAK-CENTRIC grid plots for Supplementary figure

require(ggplot2)
require(tidyverse)
require(cowplot)
require(patchwork)
theme_set(theme_cowplot())

setwd('/Users/acadiak/Desktop/CGI/Figures/Mar_2023/')

table <- read.table('/Users/acadiak/Desktop/CGI/speciesPairs/peakCentric/permutation/permutation_peakCentric_qAndFD_3x3.txt',
                    header = T) %>% tibble()

FD.table <- table %>%
  gather(key = 'combo.FD', value = 'FD', ends_with('.FD')) %>%
  separate(combo.FD, c('Peak', 'CGI', '.FD')) %>%
  select(-c(ends_with('.FD'))) %>%
  select(-c(ends_with('.q')))

q.table <- table %>%
  gather(key = 'combo.q', value = 'q', ends_with('.q')) %>%
  separate(combo.q, c('Peak', 'CGI', '.q')) %>%
  select(-c(ends_with('.q'))) %>%
  select(-c(ends_with('.FD')))

gridTable <- left_join(FD.table, q.table, 
                       by = c('speciesPair', 'speciesA', 'speciesB', 'tissue', 'mark', 'timePoint', 'Peak', 'CGI')) %>%
  mutate(timePoint = ifelse(is.na(timePoint), 'None', timePoint))

gridTable$shape <- NA
gridTable[gridTable$q<0.05&gridTable$FD>1,]$shape <- 17
gridTable[gridTable$q<0.05&gridTable$FD<1,]$shape <- 6

# store color for grid lines
gridLineColor <- 'gray80'


##############################
#### make plots for supplement

gridTable$log2FD <- log2(gridTable$FD)

function_gridPlot <- function(table, inSpeciesPair, inTissue, inMark, inTimePoint, highColor, lowLimit, highLimit) {
  x <- table %>%
    filter(speciesPair==inSpeciesPair,
           tissue==inTissue,
           mark==inMark,
           timePoint==inTimePoint) %>%
    mutate(log2FD = replace(log2FD, log2FD > highLimit, highLimit),
           log2FD = replace(log2FD, log2FD < lowLimit, lowLimit)) %>%
    ggplot(aes(y = fct_relevel(CGI, 'cgiBoth', 'cgiB', 'cgiA'),
               x = Peak,
               fill = log2FD,
               shape = shape)) +
    geom_tile(color = gridLineColor, size = 0.25) +
    theme(axis.text = element_blank(),
          axis.ticks = element_blank(),
          axis.line = element_blank(),
          legend.position = 'bottom',
          legend.title = element_blank(),
          legend.key.size = unit(0.5, 'cm')) +
    scale_fill_gradient2(low = 'gray30', high = highColor,
                         mid='white', limits=c(lowLimit, highLimit), na.value = 'red') +
    labs(x = '', y = '') +
    geom_point(aes(x = Peak, y = CGI, shape = shape), size = 1, stroke = 0.25) +
    scale_shape_identity() +
    theme(legend.position = 'none')
  return(x)
}

# Rhesus vs Mouse
pair1_brain_H3K4me3 <- function_gridPlot(gridTable, 'rheMac10_mm39', 'brain', 'H3K4me3', 'None', 'darkorange1', -2.5, 2.5)
pair1_brain_H3K27ac <- function_gridPlot(gridTable, 'rheMac10_mm39', 'brain', 'H3K27ac', 'None', 'palegreen4', -1.5, 1.5)
pair1_brain_H3K4me1 <- function_gridPlot(gridTable, 'rheMac10_mm39', 'brain', 'H3K4me1', 'None', 'mediumpurple4', -1, 1)

pair1_liver_H3K4me3 <- function_gridPlot(gridTable, 'rheMac10_mm39', 'liver', 'H3K4me3', 'None', 'darkorange1', -2.5, 2.5)
pair1_liver_H3K27ac <- function_gridPlot(gridTable, 'rheMac10_mm39', 'liver', 'H3K27ac', 'None', 'palegreen4', -1.5, 1.5)
pair1_liver_H3K4me1 <- function_gridPlot(gridTable, 'rheMac10_mm39', 'liver', 'H3K4me1', 'None', 'mediumpurple4', -1, 1)

pair1_muscle_H3K4me3 <- function_gridPlot(gridTable, 'rheMac10_mm39', 'muscle', 'H3K4me3', 'None', 'darkorange1', -2.5, 2.5)
pair1_muscle_H3K27ac <- function_gridPlot(gridTable, 'rheMac10_mm39', 'muscle', 'H3K27ac', 'None', 'palegreen4', -1.5, 1.5)
pair1_muscle_H3K4me1 <- function_gridPlot(gridTable, 'rheMac10_mm39', 'muscle', 'H3K4me1', 'None', 'mediumpurple4', -1, 1)

pair1_testis_H3K4me3 <- function_gridPlot(gridTable, 'rheMac10_mm39', 'testis', 'H3K4me3', 'None', 'darkorange1', -2.5, 2.5)
pair1_testis_H3K27ac <- function_gridPlot(gridTable, 'rheMac10_mm39', 'testis', 'H3K27ac', 'None', 'palegreen4', -1.5, 1.5)
pair1_testis_H3K4me1 <- function_gridPlot(gridTable, 'rheMac10_mm39', 'testis', 'H3K4me1', 'None', 'mediumpurple4', -1, 1)

pair1 <- (pair1_brain_H3K4me3 + pair1_brain_H3K27ac + pair1_brain_H3K4me1) / 
  (pair1_liver_H3K4me3 + pair1_liver_H3K27ac + pair1_liver_H3K4me1) /
  (pair1_muscle_H3K4me3 + pair1_muscle_H3K27ac + pair1_muscle_H3K4me1) /
  (pair1_testis_H3K4me3 + pair1_testis_H3K27ac + pair1_testis_H3K4me1) + plot_layout(nrow = 4)

ggsave('FigS26_peakCentricSupp_rhesusMouse.pdf', pair1, width = 1150, height = 1550, units = 'px')


# Rat vs Dog
pair2_brain_H3K4me3 <- function_gridPlot(gridTable, 'rn7_canFam6', 'brain', 'H3K4me3', 'None', 'darkorange1', -2.5, 2.5)
pair2_brain_H3K27ac <- function_gridPlot(gridTable, 'rn7_canFam6', 'brain', 'H3K27ac', 'None', 'palegreen4', -1.5, 1.5)
pair2_brain_H3K4me1 <- function_gridPlot(gridTable, 'rn7_canFam6', 'brain', 'H3K4me1', 'None', 'mediumpurple4', -1, 1)

pair2_liver_H3K4me3 <- function_gridPlot(gridTable, 'rn7_canFam6', 'liver', 'H3K4me3', 'None', 'darkorange1', -2.5, 2.5)
pair2_liver_H3K27ac <- function_gridPlot(gridTable, 'rn7_canFam6', 'liver', 'H3K27ac', 'None', 'palegreen4', -1.5, 1.5)
pair2_liver_H3K4me1 <- function_gridPlot(gridTable, 'rn7_canFam6', 'liver', 'H3K4me1', 'None', 'mediumpurple4', -1, 1)

pair2_muscle_H3K4me3 <- function_gridPlot(gridTable, 'rn7_canFam6', 'muscle', 'H3K4me3', 'None', 'darkorange1', -2.5, 2.5)
pair2_muscle_H3K27ac <- function_gridPlot(gridTable, 'rn7_canFam6', 'muscle', 'H3K27ac', 'None', 'palegreen4', -1.5, 1.5)
pair2_muscle_H3K4me1 <- function_gridPlot(gridTable, 'rn7_canFam6', 'muscle', 'H3K4me1', 'None', 'mediumpurple4', -1, 1)

pair2_testis_H3K4me3 <- function_gridPlot(gridTable, 'rn7_canFam6', 'testis', 'H3K4me3', 'None', 'darkorange1', -2.5, 2.5)
pair2_testis_H3K27ac <- function_gridPlot(gridTable, 'rn7_canFam6', 'testis', 'H3K27ac', 'None', 'palegreen4', -1.5, 1.5)
pair2_testis_H3K4me1 <- function_gridPlot(gridTable, 'rn7_canFam6', 'testis', 'H3K4me1', 'None', 'mediumpurple4', -1, 1)

pair2 <- (pair2_brain_H3K4me3 + pair2_brain_H3K27ac + pair2_brain_H3K4me1) / 
  (pair2_liver_H3K4me3 + pair2_liver_H3K27ac + pair2_liver_H3K4me1) /
  (pair2_muscle_H3K4me3 + pair2_muscle_H3K27ac + pair2_muscle_H3K4me1) /
  (pair2_testis_H3K4me3 + pair2_testis_H3K27ac + pair2_testis_H3K4me1) + plot_layout(nrow = 4)

ggsave('FigS26_peakCentricSupp_ratDog.pdf', pair2, width = 1150, height = 1550, units = 'px')

# Marmoset vs Cat
pair3_brain_H3K4me3 <- function_gridPlot(gridTable, 'calJac4_felCat9', 'brain', 'H3K4me3', 'None', 'darkorange1', -2.5, 2.5)
pair3_brain_H3K27ac <- function_gridPlot(gridTable, 'calJac4_felCat9', 'brain', 'H3K27ac', 'None', 'palegreen4', -1.5, 1.5)
pair3_brain_H3K4me1 <- function_gridPlot(gridTable, 'calJac4_felCat9', 'brain', 'H3K4me1', 'None', 'mediumpurple4', -1, 1)

pair3_liver_H3K4me3 <- function_gridPlot(gridTable, 'calJac4_felCat9', 'liver', 'H3K4me3', 'None', 'darkorange1', -2.5, 2.5)
pair3_liver_H3K27ac <- function_gridPlot(gridTable, 'calJac4_felCat9', 'liver', 'H3K27ac', 'None', 'palegreen4', -1.5, 1.5)
pair3_liver_H3K4me1 <- function_gridPlot(gridTable, 'calJac4_felCat9', 'liver', 'H3K4me1', 'None', 'mediumpurple4', -1, 1)

pair3_muscle_H3K4me3 <- function_gridPlot(gridTable, 'calJac4_felCat9', 'muscle', 'H3K4me3', 'None', 'darkorange1', -2.5, 2.5)
pair3_muscle_H3K27ac <- function_gridPlot(gridTable, 'calJac4_felCat9', 'muscle', 'H3K27ac', 'None', 'palegreen4', -1.5, 1.5)
pair3_muscle_H3K4me1 <- function_gridPlot(gridTable, 'calJac4_felCat9', 'muscle', 'H3K4me1', 'None', 'mediumpurple4', -1, 1)

pair3_testis_H3K4me3 <- function_gridPlot(gridTable, 'calJac4_felCat9', 'testis', 'H3K4me3', 'None', 'darkorange1', -2.5, 2.5)
pair3_testis_H3K27ac <- function_gridPlot(gridTable, 'calJac4_felCat9', 'testis', 'H3K27ac', 'None', 'palegreen4', -1.5, 1.5)
pair3_testis_H3K4me1 <- function_gridPlot(gridTable, 'calJac4_felCat9', 'testis', 'H3K4me1', 'None', 'mediumpurple4', -1, 1)

pair3 <- (pair3_brain_H3K4me3 + pair3_brain_H3K27ac + pair3_brain_H3K4me1) / 
  (pair3_liver_H3K4me3 + pair3_liver_H3K27ac + pair3_liver_H3K4me1) /
  (pair3_muscle_H3K4me3 + pair3_muscle_H3K27ac + pair3_muscle_H3K4me1) /
  (pair3_testis_H3K4me3 + pair3_testis_H3K27ac + pair3_testis_H3K4me1) + plot_layout(nrow = 4)

ggsave('FigS26_peakCentricSupp_marmosetCat.pdf', pair3, width = 1150, height = 1550, units = 'px')

# Pig vs Horse
pair4_brain_H3K4me3 <- function_gridPlot(gridTable, 'susScr11_equCab3', 'brain', 'H3K4me3', 'None', 'darkorange1', -2.5, 2.5)
pair4_brain_H3K27ac <- function_gridPlot(gridTable, 'susScr11_equCab3', 'brain', 'H3K27ac', 'None', 'palegreen4', -1.5, 1.5)
pair4_brain_H3K4me1 <- function_gridPlot(gridTable, 'susScr11_equCab3', 'brain', 'H3K4me1', 'None', 'mediumpurple4', -1, 1)

pair4_liver_H3K4me3 <- function_gridPlot(gridTable, 'susScr11_equCab3', 'liver', 'H3K4me3', 'None', 'darkorange1', -2.5, 2.5)
pair4_liver_H3K27ac <- function_gridPlot(gridTable, 'susScr11_equCab3', 'liver', 'H3K27ac', 'None', 'palegreen4', -1.5, 1.5)
pair4_liver_H3K4me1 <- function_gridPlot(gridTable, 'susScr11_equCab3', 'liver', 'H3K4me1', 'None', 'mediumpurple4', -1, 1)

pair4_muscle_H3K4me3 <- function_gridPlot(gridTable, 'susScr11_equCab3', 'muscle', 'H3K4me3', 'None', 'darkorange1', -2.5, 2.5)
pair4_muscle_H3K27ac <- function_gridPlot(gridTable, 'susScr11_equCab3', 'muscle', 'H3K27ac', 'None', 'palegreen4', -1.5, 1.5)
pair4_muscle_H3K4me1 <- function_gridPlot(gridTable, 'susScr11_equCab3', 'muscle', 'H3K4me1', 'None', 'mediumpurple4', -1, 1)

pair4_testis_H3K4me3 <- function_gridPlot(gridTable, 'susScr11_equCab3', 'testis', 'H3K4me3', 'None', 'darkorange1', -2.5, 2.5)
pair4_testis_H3K27ac <- function_gridPlot(gridTable, 'susScr11_equCab3', 'testis', 'H3K27ac', 'None', 'palegreen4', -1.5, 1.5)
pair4_testis_H3K4me1 <- function_gridPlot(gridTable, 'susScr11_equCab3', 'testis', 'H3K4me1', 'None', 'mediumpurple4', -1, 1)

pair4 <- (pair4_brain_H3K4me3 + pair4_brain_H3K27ac + pair4_brain_H3K4me1) / 
  (pair4_liver_H3K4me3 + pair4_liver_H3K27ac + pair4_liver_H3K4me1) /
  (pair4_muscle_H3K4me3 + pair4_muscle_H3K27ac + pair4_muscle_H3K4me1) /
  (pair4_testis_H3K4me3 + pair4_testis_H3K27ac + pair4_testis_H3K4me1) + plot_layout(nrow = 4)

ggsave('FigS26_peakCentricSupp_pigHorse.pdf', pair4, width = 1150, height = 1550, units = 'px')

#############################
# print results to Table S4 #
#############################

suppGridTable <- gridTable %>%
  select(speciesA, speciesB, tissue, mark, timePoint, CGI, Peak, q, log2FD) %>%
  mutate(speciesA = factor(speciesA,
                           levels = c('hg19', 'rheMac10', 'rheMac2', 'calJac4', 'mm39', 'mm9', 'rn7', 'susScr11', 'canFam6', 'felCat9'),
                           labels = c('Human', 'Rhesus', 'Rhesus', 'Marmoset', 'Mouse', 'Mouse', 'Rat', 'Pig', 'Dog', 'Cat')),
         speciesB = factor(speciesB,
                           levels = c('rheMac2', 'calJac4', 'mm39', 'mm9', 'rn7', 'susScr11', 'canFam6', 'felCat9', 'equCab3'),
                           labels = c('Rhesus', 'Marmoset', 'Mouse', 'Mouse', 'Rat', 'Pig', 'Dog', 'Cat', 'Horse')),
         tissue = factor(tissue,
                         levels = c('brain', 'liver', 'muscle', 'testis', 'devBrain', 'devLimb'),
                         labels = c('Brain', 'Liver', 'Muscle', 'Testis', 'Developing Brain', 'Developing Limb')),
         timePoint = factor(timePoint,
                            levels = c('None', '0', '1', '2', '3'),
                            labels = c('NA', '0', '1', '2', '3')),
         mark = factor(mark,
                       levels = c('H3K4me3', 'H3K27ac', 'H3K4me1', 'H3K4me2'))) %>%
  # flip labels around (Peak, CGI instead of CGI, Peak)
  pivot_wider(names_from = c(Peak, CGI),
              values_from = c(log2FD, q),
              names_sep = '.')

# write Table S4
write_delim(suppGridTable %>% arrange(mark, speciesA, speciesB, tissue),
            'TableS4_enrichment-depletion.txt',
            delim = '\t',
            quote = 'none')

# have to move some developing rows to bottom with rest, and flip order of FD and q. otherwise it's as in Table S3


