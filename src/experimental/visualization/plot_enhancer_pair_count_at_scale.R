# This program plots the distribution of number of cells the enhancer pair appears in for the
# nearly 1 million enhancer pairs that were derived from the at=scale gRNAgroup-gene pairs used in
# the Gasperini et al. paper
#
# Author: Karthik Guruvayurappan

library(ggplot2)
library(RColorBrewer)
library(rhdf5)

pair.counts <- read.csv('/iblm/netapp/data1/external/Gasperini2019/processed/at_scale_enhancer_enhancer_pairs_both_cells_count_nodups.csv')
gene.names <- h5read('/iblm/netapp/data1/external/Gasperini2019/processed/gasperini_data.h5', 'gene.names')
print(nrow(pair.counts))
pair.counts <- pair.counts[pair.counts$gene %in% gene.names, ]
pair.counts <- pair.counts[, c('enhancer_1', 'enhancer_2', 'count')]
print(nrow(pair.counts))
pair.counts$color <- 'black'
pair.counts$color[pair.counts$count > 10] <- 'red'
bins <- seq(min(pair.counts$count), max(pair.counts$count), by = 1)

plot <- ggplot(pair.counts, aes(x = count, fill = color)) +
    geom_histogram(breaks = bins) +
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
    legend.position = 'none'
    ) +
    scale_fill_manual(values = c('gray' = 'gray', 'red' = 'red'))
    
    
ggsave(
    filename = '/iblm/netapp/home/karthik/GLiMMIRS/out/23_08_18_enhancer_enhancer_at_scale_counts.pdf',
    device = 'pdf',
    plot = plot
)

ggsave(
    filename = '/iblm/netapp/home/karthik/GLiMMIRS/out/23_08_18_enhancer_enhancer_at_scale_counts.png',
    device = 'png',
    plot = plot
)

pair.counts <- pair.counts[pair.counts$count > 20, ]
print(nrow(pair.counts))

plot <- ggplot(pair.counts, aes(x = count)) +
    geom_histogram(color = 'black') +
    scale_x_continuous(expand = c(0, 0), limits = c(0, 55)) +
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
    filename = '/iblm/netapp/home/karthik/GLiMMIRS/out/23_08_16_enhancer_enhancer_at_scale_counts_20min.pdf',
    device = 'pdf',
    plot = plot
)

ggsave(
    filename = '/iblm/netapp/home/karthik/GLiMMIRS/out/23_08_16_enhancer_enhancer_at_scale_counts_20min.png',
    device = 'png',
    plot = plot
)
