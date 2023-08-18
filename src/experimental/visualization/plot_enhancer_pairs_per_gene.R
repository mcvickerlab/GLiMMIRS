# This script plots the number of enhancers pairs tested per gene for the
# enhancer pairs tested in both the smaller subset derived from Supplementary
# Table 2 of the Gasperini paper and the at-scale analysis.
#
# Author: Karthik Guruvayurappan

library(dplyr)
library(ggplot2)
library(rhdf5)

# read in enhancer pairs from supplementary table 2
enhancer.enhancer.pairs <- read.csv('/iblm/netapp/data1/external/Gasperini2019/processed/enhancer_pairs_suppl_table_2.csv')
gene.counts <- enhancer.enhancer.pairs %>% count(gene)

plot <- ggplot(gene.counts, aes(x = n)) +
    geom_histogram(color = 'black', bins = 10) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    xlab(bquote("Number of Enhancer Pairs")) + 
    ylab(bquote(Count)) +
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
    plot = plot,
    filename = '/iblm/netapp/home/karthik/GLiMMIRS/out/23_08_18_enhancer_pairs_per_gene_supp2.png',
    device = 'png'
)

# read in enhancer pairs from at-scale analysis
gene.names <- h5read('/iblm/netapp/data1/external/Gasperini2019/processed/gasperini_data.h5', 'gene.names')
enhancer.enhancer.pairs <- read.csv('/iblm/netapp/data1/external/Gasperini2019/processed/at_scale_enhancer_enhancer_pairs_both_cells_count_nodups.csv')
enhancer.enhancer.pairs <- enhancer.enhancer.pairs[enhancer.enhancer.pairs$count > 20, ]
enhancer.enhancer.pairs <- enhancer.enhancer.pairs[enhancer.enhancer.pairs$gene %in% gene.names, ]
gene.counts <- enhancer.enhancer.pairs %>% count(gene)

plot <- ggplot(gene.counts, aes(x = n)) +
    geom_histogram(color = 'black', bins = 10) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    xlab(bquote("Number of Enhancer Pairs")) + 
    ylab(bquote(Count)) +
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
    plot = plot,
    filename = '/iblm/netapp/home/karthik/GLiMMIRS/out/23_08_18_enhancer_pairs_per_gene_at_scale.png',
    device = 'png'
)
