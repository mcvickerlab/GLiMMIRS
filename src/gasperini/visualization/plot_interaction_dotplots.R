# This script plots dotplots of gene expression counts for cells with both
# enhancers perturbed for the 46 significant interactions observed in the
# analysis of 82,314 enhancer pairs.
#
# Author: Karthik Guruvayurappan

library(stats)
library(rhdf5)
library(MASS)
library(ggplot2)

# useful analysis functions
source('src/gasperini/models/analysis_helpers.R')

# create directory to hold outputs
dir.create('out/interaction_dotplots/')

# read in significant interactions and filter columns
significant_results <- get_significant_results()
significant_results <- significant_results[
  ,
  c('enhancer.1.list', 'enhancer.2.list', 'gene.list')
]
colnames(significant_results) <- c('enhancer_1', 'enhancer_2', 'gene')

# read in enhancer-guide table
enhancer_guide <- read_enhancer_guide_table()

# read in guide efficiency information
guide_info <- read_guide_efficiency_info()

# read in guide matrix
guide_matrix <- read_guide_matrix()

# read in gene expression matrix
expr_matrix <- read_expr_matrix()

# read in cell-level covariates
covariates <- read_covariates()

# iterate through significant interactions
for (i in 1:nrow(significant_results)) {

  # get name of enhancers and gene
  enhancer_1 <- significant_results[i, 'enhancer_1']
  enhancer_2 <- significant_results[i, 'enhancer_2']
  gene <- significant_results[i, 'gene']

  # get perturbation vectors corresponding to enhancers
  enh_1_perturbation <- get_enhancer_perturbation(enhancer_1, enhancer_guide,
                                                  guide_info, guide_matrix)
  enh_2_perturbation <- get_enhancer_perturbation(enhancer_2, enhancer_guide,
                                                  guide_info, guide_matrix)

  # get gene expression counts
  gene_counts <- expr_matrix[gene, ]

  # get interaction expression counts
  interaction_counts <- enh_1_perturbation * enh_2_perturbation
  interaction_counts <- gene_counts[interaction_counts > 0]
  interaction_counts <- data.frame(interaction_counts)
  interaction_counts$group <- ''

  # plot interaction dotplot
  plot <- ggplot(interaction_counts, aes(x = group, y = interaction_counts)) + 
    geom_dotplot(
      binaxis='y', 
      stackdir='center', 
      dotsize = 1.4, 
      fill = 'darkgray',
      color = 'darkgray',
    ) +
    theme_classic() +
    xlab('') +
    ylab('Expression') +
    scale_y_continuous(expand = c(0.02, 0)) +
    theme(
      axis.line = element_line(linewidth = 1),
      axis.title.y = element_text(size = 24, color = 'black'),
      axis.text = element_text(size = 24, color = 'black'),
      axis.ticks = element_line(color = 'black', linewidth = 1),
      axis.ticks.length = unit(2, 'mm'),
      plot.margin = rep(unit(10, 'mm'), 4),
      legend.text = element_text(size = 12)
    )

  ggsave(
    paste0(
        'out/interaction_dotplots/',
        enhancer_1,
        '_',
        enhancer_2,
        '_',
        gene,
        '_interaction_dotplot.pdf' 
    ),
    plot,
    device = 'pdf'
  )
}
