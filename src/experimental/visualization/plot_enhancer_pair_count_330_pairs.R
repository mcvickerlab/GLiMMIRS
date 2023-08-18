# This script plots the distribution of cells with both enhancers targeted for
# the 330 enhancer pairs derived from supplementary table 2 of the Gasperini
# paper.
#
# Author: Karthik Guruvayurappan

library(ggplot2)
library(RColorBrewer)

pair.counts <- read.csv('/iblm/netapp/data1/external/Gasperini2019/processed/enhancer_enhancer_pairs_suppl_table_2_count_both_enhancers_target.csv')
pair.counts <- pair.counts[, c('enhancer.1.list', 'enhancer.2.list', 'count.list')]

plot <- ggplot(pair.counts, aes(x = count.list)) +
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
    filename = '/iblm/netapp/home/karthik/GLiMMIRS/out/23_08_14_enhancer_enhancer_330_pair_counts.pdf',
    device = 'pdf',
    plot = plot
)

ggsave(
    filename = '/iblm/netapp/home/karthik/GLiMMIRS/out/23_08_14_enhancer_enhancer_330_pair_counts.png',
    device = 'png',
    plot = plot
)
