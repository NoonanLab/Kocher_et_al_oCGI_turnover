# 10/24/22
# Purpose: generate 3 x 3 grids for use in Fig 3, Fig 4, and Fig 6

# THIS WAS FORMERLY FIGURE 3/5/7 - NOW REVISED FIGURE 3/5/6

require(ggplot2)
require(tidyverse)
require(cowplot)
require(patchwork)
theme_set(theme_cowplot())

setwd('/Users/acadiak/Desktop/CGI/Figures/Mar_2023/')

table <- read.table('/Users/acadiak/Desktop/CGI/speciesPairs/CGIcentric/permutation/permutation_qAndFD_3x3.txt',
                    header = T) %>% tibble()

FD.table <- table %>%
  gather(key = 'combo.FD', value = 'FD', ends_with('.FD')) %>%
  separate(combo.FD, c('CGI', 'Peak', '.FD')) %>%
  select(-c(ends_with('.FD'))) %>%
  select(-c(ends_with('.q')))

q.table <- table %>%
  gather(key = 'combo.q', value = 'q', ends_with('.q')) %>%
  separate(combo.q, c('CGI', 'Peak', '.q')) %>%
  select(-c(ends_with('.q'))) %>%
  select(-c(ends_with('.FD')))

gridTable <- left_join(FD.table, q.table, 
                       by = c('speciesPair', 'speciesA', 'speciesB', 'tissue', 'mark', 'timePoint', 'CGI', 'Peak')) %>%
  mutate(timePoint = ifelse(is.na(timePoint), 'None', timePoint))

gridTable$shape <- NA
gridTable[gridTable$q<0.05&gridTable$FD>1,]$shape <- 17
gridTable[gridTable$q<0.05&gridTable$FD<1,]$shape <- 6

# add column for log2(FD)
gridTable$log2FD <- log2(gridTable$FD)

# get example values for text
gridTable %>% filter(speciesPair == 'rheMac10_mm39', tissue == 'brain', mark == 'H3K4me3')

# store color for grid lines
gridLineColor <- 'gray80'

# make grids for Fig 3C
Fig3_H3K4me3 <- gridTable %>%
  filter(speciesPair=='rheMac10_mm39',
         tissue=='brain',
         mark=='H3K4me3') %>%
  ggplot(aes(y = fct_relevel(Peak, 'peakBoth', 'peakB', 'peakA'),
             x = CGI,
             fill = log2(FD),
             shape = shape)) +
  geom_tile(color = gridLineColor, size = 0.5) +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        legend.position = 'bottom',
        legend.title = element_blank(),
        legend.key.size = unit(0.5, 'cm')) +
  scale_fill_gradient2(low = 'gray30', high = 'darkorange', 
                       mid='white', limits=c(-2.5,2.5),na.value='red') +
  labs(x = '', y = '') +
  geom_point(aes(x = CGI, y = Peak, shape = shape), size = 2) +
  scale_shape_identity()

Fig3_H3K27ac <- gridTable %>%
  filter(speciesPair=='rheMac10_mm39',
         tissue=='brain',
         mark=='H3K27ac') %>%
  ggplot(aes(y = fct_relevel(Peak, 'peakBoth', 'peakB', 'peakA'),
             x = CGI,
             fill = log2(FD),
             shape = shape)) +
  geom_tile(color = gridLineColor, size = 0.5) +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        legend.position = 'bottom',
        legend.title = element_blank(),
        legend.key.size = unit(0.5, 'cm')) +
  scale_fill_gradient2(low = 'gray30', high = 'palegreen4', 
                       mid='white', limits=c(-1.5,1.5),na.value='red') +
  labs(x = '', y = '') +
  geom_point(aes(x = CGI, y = Peak, shape = shape), size = 2) +
  scale_shape_identity()

Fig3_H3K4me1 <- gridTable %>%
  filter(speciesPair=='rheMac10_mm39',
         tissue=='brain',
         mark=='H3K4me1') %>%
  ggplot(aes(y = fct_relevel(Peak, 'peakBoth', 'peakB', 'peakA'),
             x = CGI,
             fill = log2(FD),
             shape = shape)) +
  geom_tile(color = gridLineColor, size = 0.5) +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        legend.position = 'bottom',
        legend.title = element_blank(),
        legend.key.size = unit(0.5, 'cm')) +
  scale_fill_gradient2(low = 'gray30', high = 'mediumpurple4', 
                       mid='white', limits=c(-1,1),na.value='red') +
  labs(x = '', y = '') +
  geom_point(aes(x = CGI, y = Peak, shape = shape), size = 2) +
  scale_shape_identity()


Fig3D <- Fig3_H3K4me3 + Fig3_H3K27ac + Fig3_H3K4me1 + plot_layout(ncol = 3)

ggsave('Fig3C_exampleGridplots.pdf', Fig3D, width = 1825, height = 850, units = 'px')


###### grids for Fig 3D
function_gridPlot <- function(table, inSpeciesPair, inTissue, inMark, inTimePoint, highColor, limit) {
  x <- table %>%
    filter(speciesPair==inSpeciesPair,
           tissue==inTissue,
           mark==inMark,
           timePoint==inTimePoint) %>%
    ggplot(aes(y = fct_relevel(Peak, 'peakBoth', 'peakB', 'peakA'),
               x = CGI,
               fill = log2(FD),
               shape = shape)) +
    geom_tile(color = gridLineColor, size = 0.25) +
    theme(axis.text = element_blank(),
          axis.ticks = element_blank(),
          axis.line = element_blank(),
          legend.position = 'bottom',
          legend.title = element_blank(),
          legend.key.size = unit(0.5, 'cm')) +
    scale_fill_gradient2(low = 'gray30', high = highColor, 
                         mid='white', limits=c(-limit, limit),na.value='red') +
    labs(x = '', y = '') +
    geom_point(aes(x = CGI, y = Peak, shape = shape), size = 1, stroke = 0.25) +
    scale_shape_identity() +
    theme(legend.position = 'none')
  return(x)
}


Fig3D_pair1_liver_H3K4me3 <- function_gridPlot(gridTable, 'rn7_canFam6', 'liver', 'H3K4me3', 'None', 'darkorange1', 2.5)
Fig3D_pair1_liver_H3K27ac <- function_gridPlot(gridTable, 'rn7_canFam6', 'liver', 'H3K27ac', 'None', 'palegreen4', 1.5)
Fig3D_pair1_liver_H3K4me1 <- function_gridPlot(gridTable, 'rn7_canFam6', 'liver', 'H3K4me1', 'None', 'mediumpurple4', 1)

Fig3D_pair1_muscle_H3K4me3 <- function_gridPlot(gridTable, 'rn7_canFam6', 'muscle', 'H3K4me3', 'None', 'darkorange1', 2.5)
Fig3D_pair1_muscle_H3K27ac <- function_gridPlot(gridTable, 'rn7_canFam6', 'muscle', 'H3K27ac', 'None', 'palegreen4', 1.5)
Fig3D_pair1_muscle_H3K4me1 <- function_gridPlot(gridTable, 'rn7_canFam6', 'muscle', 'H3K4me1', 'None', 'mediumpurple4', 1)

Fig3D_pair1 <- (Fig3D_pair1_liver_H3K4me3 + Fig3D_pair1_liver_H3K27ac + Fig3D_pair1_liver_H3K4me1) / (Fig3D_pair1_muscle_H3K4me3 + Fig3D_pair1_muscle_H3K27ac + Fig3D_pair1_muscle_H3K4me1) +
  plot_layout(nrow = 2)

ggsave('Fig3D_pair1_ratDog.pdf', Fig3D_pair1, width = 1150, height = 800, units = 'px')


Fig4A_pair1_devbrain_H3K27ac <- function_gridPlot(gridTable, 'hg19_mm9', 'devBrain', 'H3K27ac', 1, 'palegreen4', 1.5)
Fig4A_pair1_devbrain_H3K4me2 <- function_gridPlot(gridTable, 'hg19_mm9', 'devBrain', 'H3K4me2', 1, 'darkgoldenrod1', 2)
Fig4A_pair1_devlimb_H3K27ac <- function_gridPlot(gridTable, 'hg19_mm9', 'devLimb', 'H3K27ac', 1, 'palegreen4', 1.5)

# repeat bottom - remove in illustrator
Fig4A_pair1 <- (Fig4A_pair1_devbrain_H3K27ac + Fig4A_pair1_devbrain_H3K4me2) / (Fig4A_pair1_devlimb_H3K27ac + Fig4A_pair1_devlimb_H3K27ac)

ggsave('Fig4A_pair1_humanMouse.pdf', Fig4A_pair1, width = 790, height = 800, units = 'px')

############################
######## Supplement ########
############################

# add column for log2(FD)
gridTable$log2FD <- log2(gridTable$FD)

# store function for making triangle grid plots of all species pairs
# speciesA on y, speciesB on x
function_allSpecies_BvsA <- function(table, inTissue, inMark, inTimePoint, highColor, lowLimit, highLimit) {
  x <- table %>%
  filter(tissue==inTissue,
         mark==inMark,
         timePoint == inTimePoint) %>%
  mutate(log2FD = replace(log2FD, log2FD > highLimit, highLimit),
         log2FD = replace(log2FD, log2FD < lowLimit, lowLimit)) %>%
  ggplot(aes(y = fct_relevel(Peak, 'peakBoth', 'peakB', 'peakA'),
             x = CGI,
             fill = log2FD,
             shape = shape)) +
  facet_grid(speciesB ~ speciesA, switch = 'both') +
  geom_tile(color = gridLineColor, size = 0.25) +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        legend.position = 'bottom',
        legend.title = element_blank(),
        legend.key.size = unit(0.5, 'cm'),
        panel.spacing = unit(0.25, 'lines'),
        strip.background = element_blank(),
        strip.text = element_text(size = 10),
        strip.text.y.left = element_text(angle = 0),
        plot.title = element_text(hjust = 0.3)) +
  scale_fill_gradient2(low = 'gray30', high = highColor, 
                       mid='white', limits=c(lowLimit, highLimit), na.value='red') +
  labs(x = 'Species A', y = 'Species B', title = paste0(toupper(substr(inTissue, 0, 1)), substr(inTissue, 2, 8))) +
  geom_point(aes(x = CGI, y = Peak, shape = shape), size = 1, stroke = 0.25) +
  scale_shape_identity()
  return(x)
}

cp.marks.alternateShades <- c('darkorange1', 'tan1', 
                              'palegreen4', 'darkseagreen3',
                              'mediumpurple4', 'mediumpurple2',
                              'darkgoldenrod1', 'gold')

# store Roller data separately and re-factor genome labels into species names
Roller_gridTable <- gridTable %>%
  filter(speciesA %in% c('rheMac10', 'calJac4', 'mm39', 'rn7', 'susScr11', 'canFam6', 'felCat9'),
         mark %in% c('H3K4me3', 'H3K27ac', 'H3K4me1')) %>%
  mutate(speciesA = factor(speciesA,
                           levels = c('rheMac10', 'calJac4', 'mm39', 'rn7', 'susScr11', 'canFam6', 'felCat9'),
                           labels = c('Rhesus', 'Marmoset', 'Mouse', 'Rat', 'Pig', 'Dog', 'Cat')),
         speciesB = factor(speciesB,
                           levels = c('calJac4', 'mm39', 'rn7', 'susScr11', 'canFam6', 'felCat9', 'equCab3'),
                           labels = c('Marmoset', 'Mouse', 'Rat', 'Pig', 'Dog', 'Cat', 'Horse')))

# make plots
H3K4me3_brain <- function_allSpecies_BvsA(Roller_gridTable, 'brain', 'H3K4me3', 'None', 'darkorange1', -2.5, 2.5) # 3 down
H3K4me3_liver <- function_allSpecies_BvsA(Roller_gridTable, 'liver', 'H3K4me3', 'None', 'darkorange1', -2.5, 2.5) # ~10 down
H3K4me3_muscle <- function_allSpecies_BvsA(Roller_gridTable, 'muscle', 'H3K4me3', 'None', 'darkorange1', -2.5, 2.5) # 1 down
H3K4me3_testis <- function_allSpecies_BvsA(Roller_gridTable, 'testis', 'H3K4me3', 'None', 'darkorange1', -2.5, 2.5) # 4 down

H3K27ac_brain <- function_allSpecies_BvsA(Roller_gridTable, 'brain', 'H3K27ac', 'None', 'palegreen4', -1.5, 1.5)
H3K27ac_liver <- function_allSpecies_BvsA(Roller_gridTable, 'liver', 'H3K27ac', 'None', 'palegreen4', -1.5, 1.5)
H3K27ac_muscle <- function_allSpecies_BvsA(Roller_gridTable, 'muscle', 'H3K27ac', 'None', 'palegreen4', -1.5, 1.5)
H3K27ac_testis <- function_allSpecies_BvsA(Roller_gridTable, 'testis', 'H3K27ac', 'None', 'palegreen4', -1.5, 1.5)

H3K4me1_brain <- function_allSpecies_BvsA(Roller_gridTable, 'brain', 'H3K4me1', 'None', 'mediumpurple4', -1, 1) # 1 down
H3K4me1_liver <- function_allSpecies_BvsA(Roller_gridTable, 'liver', 'H3K4me1', 'None', 'mediumpurple4', -1, 1)
H3K4me1_muscle <- function_allSpecies_BvsA(Roller_gridTable, 'muscle', 'H3K4me1', 'None', 'mediumpurple4', -1, 1)
H3K4me1_testis <- function_allSpecies_BvsA(Roller_gridTable, 'testis', 'H3K4me1', 'None', 'mediumpurple4', -1, 1) # ~12 down, 4 up

# now Fig S23-24
H3K4me3 <- H3K4me3_brain + H3K4me3_muscle + H3K4me3_liver + H3K4me3_testis + plot_layout(ncol = 2)
ggsave('FigS18-19_H3K4me3_allCombos.pdf', H3K4me3, height = 3500, width = 3300, units = 'px')

# now Fig S25-26
H3K27ac <- H3K27ac_brain + H3K27ac_muscle + H3K27ac_liver + H3K27ac_testis + plot_layout(ncol = 2)
ggsave('FigS20-21_H3K27ac_allCombos.pdf', H3K27ac, height = 3500, width = 3300, units = 'px')

# now Fig S27-28
H3K4me1 <- H3K4me1_brain + H3K4me1_muscle + H3K4me1_liver + H3K4me1_testis + plot_layout(ncol = 2)
ggsave('FigS22-23_H3K4me1_allCombos.pdf', H3K4me1, height = 3500, width = 3300, units = 'px')


# store Noonan data separately and re-factor genome labels into species names
Noonan_gridTable <- gridTable %>%
  filter(speciesA %in% c('hg19', 'rheMac2')) %>%
  mutate(speciesA = factor(speciesA,
                           levels = c('hg19', 'rheMac2'),
                           labels = c('Human', 'Rhesus')),
         speciesB = factor(speciesB,
                           levels = c('rheMac2', 'mm9'),
                           labels = c('Rhesus', 'Mouse')))

# make plots for Noonan
H3K27ac_devBrain_0 <- function_allSpecies_BvsA(Noonan_gridTable, 'devBrain', 'H3K27ac', '0', 'palegreen4', -1.5, 1.5)
H3K27ac_devBrain_1 <- function_allSpecies_BvsA(Noonan_gridTable, 'devBrain', 'H3K27ac', '1', 'palegreen4', -1.5, 1.5)
H3K27ac_devBrain_2 <- function_allSpecies_BvsA(Noonan_gridTable, 'devBrain', 'H3K27ac', '2', 'palegreen4', -1.5, 1.5)
H3K27ac_devBrain_3 <- function_allSpecies_BvsA(Noonan_gridTable, 'devBrain', 'H3K27ac', '3', 'palegreen4', -1.5, 1.5)

H3K27ac_devLimb_0 <- function_allSpecies_BvsA(Noonan_gridTable, 'devLimb', 'H3K27ac', '0', 'palegreen4', -1.5, 1.5)
H3K27ac_devLimb_1 <- function_allSpecies_BvsA(Noonan_gridTable, 'devLimb', 'H3K27ac', '1', 'palegreen4', -1.5, 1.5)
H3K27ac_devLimb_2 <- function_allSpecies_BvsA(Noonan_gridTable, 'devLimb', 'H3K27ac', '2', 'palegreen4', -1.5, 1.5)
H3K27ac_devLimb_3 <- function_allSpecies_BvsA(Noonan_gridTable, 'devLimb', 'H3K27ac', '3', 'palegreen4', -1.5, 1.5)

H3K4me2_devBrain_0 <- function_allSpecies_BvsA(Noonan_gridTable, 'devBrain', 'H3K4me2', '0', 'darkgoldenrod1', -2, 2)
H3K4me2_devBrain_1 <- function_allSpecies_BvsA(Noonan_gridTable, 'devBrain', 'H3K4me2', '1', 'darkgoldenrod1', -2, 2)
H3K4me2_devBrain_2 <- function_allSpecies_BvsA(Noonan_gridTable, 'devBrain', 'H3K4me2', '2', 'darkgoldenrod1', -2, 2)
H3K4me2_devBrain_3 <- function_allSpecies_BvsA(Noonan_gridTable, 'devBrain', 'H3K4me2', '3', 'darkgoldenrod1', -2, 2)

# print
devBrain_H3K27ac <- H3K27ac_devBrain_0 + H3K27ac_devBrain_1 + H3K27ac_devBrain_2 + H3K27ac_devBrain_3 + plot_layout(ncol = 4)
devBrain_H3K4me2 <- H3K4me2_devBrain_0 + H3K4me2_devBrain_1 + H3K4me2_devBrain_2 + H3K4me2_devBrain_3 + plot_layout(ncol = 4)
devLimb_H3K27ac <- H3K27ac_devLimb_0 + H3K27ac_devLimb_1 + H3K27ac_devLimb_2 + H3K27ac_devLimb_3 + plot_layout(ncol = 4)

# now Fig S33
ggsave('FigS27_Noonan_devBrain_H3K27ac.pdf', devBrain_H3K27ac, height = 900, width = 2700, units = 'px')
ggsave('FigS27_Noonan_devBrain_H3K4me2.pdf', devBrain_H3K4me2, height = 900, width = 2700, units = 'px')
ggsave('FigS27_Noonan_devLimb_H3K27ac.pdf', devLimb_H3K27ac, height = 900, width = 2700, units = 'px')


#################
# LIVER TF DATA #
#################

# add column for log2(FD)
gridTable$log2FD <- log2(gridTable$FD)

function_allSpecies_BvsA_TF <- function(table, inTissue, inMark, inTimePoint, highColor, lowLimit, highLimit) {
  x <- table %>%
    filter(tissue==inTissue,
           mark==inMark,
           timePoint == inTimePoint) %>%
    mutate(log2FD = replace(log2FD, log2FD > highLimit, highLimit),
           log2FD = replace(log2FD, log2FD < lowLimit, lowLimit)) %>%
    ggplot(aes(y = fct_relevel(Peak, 'peakBoth', 'peakB', 'peakA'),
               x = CGI,
               fill = log2FD,
               shape = shape)) +
    facet_grid(speciesB ~ speciesA, switch = 'both') +
    geom_tile(color = gridLineColor, size = 0.25) +
    theme(axis.text = element_blank(),
          axis.ticks = element_blank(),
          axis.line = element_blank(),
          legend.position = 'bottom',
          legend.title = element_blank(),
          legend.key.size = unit(0.5, 'cm'),
          panel.spacing = unit(0.25, 'lines'),
          strip.background = element_blank(),
          strip.text = element_text(size = 10),
          strip.text.y.left = element_text(angle = 0),
          plot.title = element_text(hjust = 0.3)) +
    scale_fill_gradient2(low = 'gray30', high = highColor, 
                         mid='white', limits=c(lowLimit, highLimit), na.value='red') +
    labs(x = 'Species A', y = 'Species B', title = inMark) +
    geom_point(aes(x = CGI, y = Peak, shape = shape), size = 1, stroke = 0.25) +
    scale_shape_identity()
  return(x)
}

LiverTF_gridTable <- gridTable %>%
  filter(mark %in% c('CTCF', 'CEBPA', 'FOXA1', 'HNF4A', 'HNF6')) %>%
  mutate(speciesA = factor(speciesA,
                           levels = c('rheMac10', 'mm39', 'rn7'),
                           labels = c('Rhesus', 'Mouse', 'Rat')),
         speciesB = factor(speciesB,
                           levels = c('mm39', 'rn7', 'canFam6'),
                           labels = c('Mouse', 'Rat', 'Dog')))

# set functions that automatically do max and min values
CTCF <- function_allSpecies_BvsA_TF(LiverTF_gridTable, 'liver', 'CTCF', 'None', 'turquoise4',
                                    floor(min(LiverTF_gridTable %>% filter(mark == 'CTCF') %>% select(log2FD)) / 0.5) * 0.5,
                                    ceiling(max(LiverTF_gridTable %>% filter(mark == 'CTCF') %>% select(log2FD)) / 0.5) * 0.5)
CEBPA <- function_allSpecies_BvsA_TF(LiverTF_gridTable, 'liver', 'CEBPA', 'None', 'seagreen',
                                     floor(min(LiverTF_gridTable %>% filter(mark == 'CEBPA') %>% select(log2FD)) / 0.5) * 0.5,
                                     ceiling(max(LiverTF_gridTable %>% filter(mark == 'CEBPA') %>% select(log2FD)) / 0.5) * 0.5)
FOXA1 <- function_allSpecies_BvsA_TF(LiverTF_gridTable, 'liver', 'FOXA1', 'None', 'firebrick3',
                                     floor(min(LiverTF_gridTable %>% filter(mark == 'FOXA1') %>% select(log2FD)) / 0.5) * 0.5,
                                     ceiling(max(LiverTF_gridTable %>% filter(mark == 'FOXA1') %>% select(log2FD)) / 0.5) * 0.5)
HNF4A <- function_allSpecies_BvsA_TF(LiverTF_gridTable, 'liver', 'HNF4A', 'None', 'magenta4',
                                     floor(min(LiverTF_gridTable %>% filter(mark == 'HNF4A') %>% select(log2FD)) / 0.5) * 0.5,
                                     ceiling(max(LiverTF_gridTable %>% filter(mark == 'HNF4A') %>% select(log2FD)) / 0.5) * 0.5)
HNF6 <- function_allSpecies_BvsA_TF(LiverTF_gridTable, 'liver', 'HNF6', 'None', 'royalblue4',
                                    floor(min(LiverTF_gridTable %>% filter(mark == 'HNF6') %>% select(log2FD)) / 0.5) * 0.5,
                                    ceiling(max(LiverTF_gridTable %>% filter(mark == 'HNF6') %>% select(log2FD)) / 0.5) * 0.5)

# now Fig S46
AllTFs <- CTCF + HNF6 + HNF4A + FOXA1 + CEBPA + plot_layout(nrow = 5)
ggsave('FigS39_allTFs.pdf', AllTFs, height = 4800, width = 800, units = 'px')

#############################
# print results to Table S3 #
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
                       levels = c('H3K4me3', 'H3K27ac', 'H3K4me1', 'H3K4me2', 'CTCF', 'HNF4A', 'HNF6', 'FOXA1', 'CEBPA'))) %>%
  pivot_wider(names_from = c(CGI, Peak),
              values_from = c(log2FD, q),
              names_sep = '.')

# write Table S3
write_delim(suppGridTable %>% arrange(mark, speciesA, speciesB, tissue),
            'TableS3_enrichment-depletion.txt',
            delim = '\t',
            quote = 'none')

# have to move some developing rows to bottom with rest. otherwise it's as in Table S3


