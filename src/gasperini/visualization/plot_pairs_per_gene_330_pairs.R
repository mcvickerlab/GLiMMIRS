# This script plots the distribution of number of enhancer pairs per gene for
# the set of 264 tested enhancer pairs.
#
# Author: Karthik Guruvayurappan

library(stats)
library(dplyr)
library(ggplot2)

# read in results from 330 enhancer-enhancer-gene sets
results <- read.csv('data/gasperini/processed/enhancer_pairs_330_models.csv')

# filter out NA values
results <- results[complete.cases(results), ]

# plot distribution of number of enhancer pairs per gene
plot_df <- results %>% count(gene.list)

plot <- ggplot(plot_df, aes(x = n)) +
  geom_histogram(fill = 'darkgray', color = 'darkgray') +
  theme_classic() +
  scale_x_continuous(expand = c(0.02, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  xlab("Number of Enhancer Pairs") + 
  ylab("Count") +
  theme(
    axis.line = element_line(linewidth = 1),
    axis.title.x = element_text(size = 24, color = 'black'),
    axis.title.y = element_text(size = 24, color = 'black'),
    axis.text = element_text(size = 24, color = 'black'),
    axis.ticks = element_line(color = 'black', linewidth = 1),
    axis.ticks.length = unit(2, 'mm'),
    plot.margin = rep(unit(10, 'mm'), 4),
  )

ggsave(
  plot = plot,
  filename = 'out/pairs_per_gene_330_pairs.pdf',
  device = 'pdf'
)

