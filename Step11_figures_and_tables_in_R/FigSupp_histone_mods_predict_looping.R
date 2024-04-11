# 3/20/24
# Purpose: see if histone mods are predictive of looping changes using published liver data

# MAKES REVISED FIG S4

require(ggplot2)
require(cowplot)
require(tidyverse)
theme_set(theme_cowplot())

setwd('/Users/acadiak/Desktop/Yale/!Writing/Kocher_CpG/GB_Revisions/Revised_code/Step12_revisions/Loops')

filePathToData <- '/Users/acadiak/Desktop/Yale/!Writing/Kocher_CpG/GB_Revisions/Revised_code/Step12_revisions/Loops/summaryTables_pad1kb'

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
  separate(filename, into = c('mark'),
           remove = T, sep = '_', extra = 'drop')

# store order for tissues and marks and factor in that order
marks <- c('H3K27ac','H3K4me3','H3K4me2','H3K4me1')

# store binary version of whether there is an interaction
#table$interactionBinary <- 0
#table[table$numberInteractions>0,]$interactionBinary <- 1
  
# refactor
table$mark <- factor(table$mark, levels = marks)
table$CGI <- factor(table$CGI, levels = c(0, 1), labels = c('Without oCGI', 'With oCGI'))
table$peakBinary <- factor(table$peakBinary, levels = c(1, 0), labels = c('+', '-'))

# do wilcoxon test
wilcox <- table %>%
  filter(grepl('dELS', cCRE)) %>%
  group_by(mark, CGI) %>%
  do(w = wilcox.test(numberInteractions~peakBinary, data=., paired = F)) %>%
  summarize(mark, CGI, Wilcox = w$p.value)

# adjust
wilcox$Wilcox.adj <- p.adjust(wilcox$Wilcox, method = 'BH')

# add columns for plotting with geom_segment and geom_text
wilcox$x <- 1.1
wilcox$xend <- 1.9
wilcox$y <- 7
wilcox$yend <- wilcox$y
wilcox$stars <- ''
wilcox[wilcox$Wilcox.adj<0.05,]$stars <- '*'
wilcox$stars.x <- (wilcox$x + wilcox$xend) / 2

# drop rows without a star to prevent line from being plotted (it's all of them)
wilcox <- wilcox[wilcox$Wilcox.adj < 0.05,]

# plot with colors
table$color <- paste0(table$mark, " ", table$CGI, " ", table$peakBinary)
table$color <- factor(table$color,
                      levels = c('H3K27ac With oCGI -', 'H3K27ac With oCGI +', 'H3K27ac Without oCGI -', 'H3K27ac Without oCGI +',
                                 'H3K4me3 With oCGI -', 'H3K4me3 With oCGI +', 'H3K4me3 Without oCGI -', 'H3K4me3 Without oCGI +',
                                 'H3K4me2 With oCGI -', 'H3K4me2 With oCGI +', 'H3K4me2 Without oCGI -', 'H3K4me2 Without oCGI +',
                                 'H3K4me1 With oCGI -', 'H3K4me1 With oCGI +', 'H3K4me1 Without oCGI -', 'H3K4me1 Without oCGI +'))
cp.marks.alternateShades <- c('gray80', 'palegreen4', 'gray80', 'darkseagreen3',
                              'gray80', 'darkorange1','gray80', 'tan1', 
                              'gray80', 'darkgoldenrod1','gray80', 'gold',
                              'gray80', 'mediumpurple4','gray80', 'mediumpurple2')

# define function for getting 90% confidence interval, 25% quantile, median, 75% quantile
quantiles_90 <- function(x) {
  r <- quantile(x, probs=c(0.10, 0.25, 0.5, 0.75, 0.90))
  names(r) <- c("ymin", "lower", "middle", "upper", "ymax")
  r
}

# boxplots
boxplots <- table %>%
  filter(grepl('dELS', cCRE)) %>%
  ggplot(aes(x = peakBinary, y = numberInteractions, fill = color)) + 
  facet_grid(mark ~ CGI, scales = 'free') +
  stat_summary(fun.data = quantiles_90, 
               geom='boxplot', 
               position=position_dodge2(preserve = 'single'),
               lwd = 0.25, width = 0.75) +
  scale_fill_manual(values = cp.marks.alternateShades) +
  coord_cartesian(ylim = c(0, 11)) +
  theme_bw() + 
  theme(panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12),
        axis.title = element_blank(),
        legend.position = 'none') +
  labs(y = 'Number of Interactions',
       x = 'Peak') +
  geom_segment(data = wilcox, 
               aes(x = x, xend = xend, 
                   y = y, yend = yend), 
               inherit.aes = F,
               lwd = 0.5) +
  geom_text(data = wilcox,
            aes(x = stars.x , y = y + 0.3, label = stars), inherit.aes = F, size = 6)

ggsave('FigS4_marks_predict_loops.pdf', boxplots, height = 1800, width = 1300, units = 'px', limitsize = F)

