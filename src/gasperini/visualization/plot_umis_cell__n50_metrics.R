# This script plots a histogram of the log10(nUMIs/cells with a UMI) per gene
# to determine which genes have a high expression count in only a few cells.
#
# Author: Karthik Guruvayurappan

library(ggplot2)
source('src/experimental/visualization/plotting_helpers.R')
source('src/experimental/models/analysis_helpers.R')

# read in metric
umi_cells_metric <- read.csv(
  'data/experimental/processed/umi_cells_n50_metrics.csv'
)
colnames(umi_cells_metric) <- c('gene', 'cell_umi_metric', 'n50')

# read in significant results
significant_results <- get_significant_results()

# get genes with significant interaction
significant_genes <- unique(significant_results$gene.list)

# add colors to data frame
umi_cells_metric$significant_gene <- FALSE
umi_cells_metric$significant_gene[
  umi_cells_metric$gene %in% significant_genes
] <- TRUE

# log-transform N50-style metric
umi_cells_metric$n50 <- log2(umi_cells_metric$n50 + 1)

# plot scatterplot of the two metrics
plot <- plot_scatterplot(umi_cells_metric, cell_umi_metric, n50) +
  geom_point(aes(color = significant_gene, alpha = 0.1)) +
  scale_color_manual(
    values = c('black', 'red'),
    breaks = c(FALSE, TRUE)
  )

ggsave(
  filename = 'out/umi_cell_n50_scatterplot.png',
  device = 'png',
  plot = plot
)

# read in expression matrix
h5_name <- 'data/experimental/processed/gasperini_data.h5'
expr_matrix <- read_expr_matrix(h5_name)

# count UMIs per gene
umis_per_gene <- rowSums(expr_matrix)

# add to plotting dataframe
umi_cells_metric$numis <- umis_per_gene
umi_cells_metric$numis <- log10(umi_cells_metric$numis)

plot <- plot_scatterplot(umi_cells_metric, numis, n50) +
  geom_point(
    data = subset(umi_cells_metric, significant_gene == FALSE), 
    aes(color = significant_gene)
  ) +
  geom_point(
    data = subset(umi_cells_metric, significant_gene == TRUE), 
    aes(color = significant_gene)
  ) + 
  scale_color_manual(
    values = c('black', 'red'),
    breaks = c(FALSE, TRUE)
  )

ggsave(
  filename = 'out/gene_umis_n50_scatterplot.png',
  device = 'png',
  plot = plot
)

plot <- plot_scatterplot(umi_cells_metric, numis, cell_umi_metric) +
  geom_point(
    data = subset(umi_cells_metric, significant_gene == FALSE), 
    aes(color = significant_gene)
  ) +
  geom_point(
    data = subset(umi_cells_metric, significant_gene == TRUE), 
    aes(color = significant_gene)
  ) + 
  scale_color_manual(
    values = c('black', 'red'),
    breaks = c(FALSE, TRUE)
  )

ggsave(
  filename = 'out/gene_umis_cell_metric_scatterplot.png',
  device = 'png',
  plot = plot
)


