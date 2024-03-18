# This script plots the distribution of number of cells with both enhancers
# perturbed for the enhancer pairs derived by using enhancer-gene links below
# an FDR threshold of 20%.
#
# Author: Karthik Guruvayurappan

library(ggplot2)
library(RColorBrewer)
library(rhdf5)

pair.counts <- read.csv('data/experimental/processed/enhancer_pair_perturbation_counts_FDR20.csv')
print(nrow(pair.counts))
pair.counts$color <- 'black'
pair.counts$color[pair.counts$count.list > 10] <- 'red'
bins <- seq(min(pair.counts$count.list), max(pair.counts$count.list), by = 1)

plot <- ggplot(pair.counts, aes(x = count.list, fill = color)) +
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
    filename = '/iblm/netapp/home/karthik/GLiMMIRS/out/23_08_18_enhancer_enhancer_FDR20_counts.pdf',
    device = 'pdf',
    plot = plot
)

ggsave(
    filename = '/iblm/netapp/home/karthik/GLiMMIRS/out/23_08_18_enhancer_enhancer_FDR20_counts.png',
    device = 'png',
    plot = plot
)
