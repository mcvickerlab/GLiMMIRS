# This script plots a histogram of the log10(nUMIs/cells with a UMI) per gene
# to determine which genes have a high expression count in only a few cells.
#
# Author: Karthik Guruvayurappan

library(ggplot2)
source('src/experimental/visualization/plotting_helpers.R')
source('src/experimental/models/analysis_helpers.R')

# read in metric
umi_cells_metric <- read.csv(
  'data/experimental/processed/umi_cells_metric.csv'
)
colnames(umi_cells_metric) <- c('gene', 'metric')

plot <- plot_histogram(umi_cells_metric, metric)

ggsave(
  filename = 'out/umi_cells_metric.png',
  device = 'png',
  plot = plot
)

# get significant results and genes
significant_results <- get_significant_results()

# get significant genes
significant_genes <- unique(significant_results$gene.list)

# get UMI cells metric for significant genes
umi_cells_metric_significant_genes <- umi_cells_metric[
  umi_cells_metric$gene %in% significant_genes,

]

# plot a histogram for the significant genes
plot <- plot_histogram(umi_cells_metric_significant_genes, metric)

ggsave(
  filename = 'out/umi_cells_metric_interaction_genes.png',
  device = 'png',
  plot = plot
)





