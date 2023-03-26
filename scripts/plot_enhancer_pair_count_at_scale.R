# This program plots the distribution of number of cells the enhancer pair appears in for the
# nearly 1 million enhancer pairs that were derived from the at=scale gRNAgroup-gene pairs used in
# the Gasperini et al. paper
#
# Author: Karthik Guruvayurappan

library(ggplot2)
library(RColorBrewer)

pair.counts <- read.csv('/iblm/netapp/data1/external/Gasperini2019/processed/at_scale_enhancer_enhancer_pairs_both_cells_count_nodups.csv')
pair.counts <- pair.counts[, c('enhancer_1', 'enhancer_2', 'count')]

plot <- ggplot(pair.counts, aes(x = count)) +
    geom_histogram(color = 'black') +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    xlab(bquote("Number of Cells")) + 
    ylab(bquote("Count of Enhancer Pairs")) +
    theme_classic() +
    theme(
    axis.line = element_line(linewidth = 1),
    axis.title.x = element_text(size = 20, color = 'black'),
    axis.title.y = element_text(size = 20, color = 'black'),
    axis.text = element_text(size = 20, color = 'black'),
    axis.ticks = element_line(color = 'black', linewidth = 1),
    axis.ticks.length = unit(2, 'mm'),
    plot.margin = rep(unit(10, 'mm'), 4),
    )
    
ggsave(
    filename = '/iblm/netapp/home/karthik/GLiMMIRS/plots/23_03_25_enhancer_enhancer_at_scale_counts.pdf',
    device = 'pdf',
    plot = plot
)
