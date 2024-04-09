# This script plots a boxplot of the MOI of cells with the double perturbation
# compared to all other cells for the significant interactions observed in the
# Gasperini dataset.
#
# Author: Karthik Guruvayurappan

library(rhdf5)
library(MASS)
library(stats)
library(ggplot2)

# source helper functions for analysis
source('src/gasperini/models/analysis_helpers.R')

# create output directory
dir.create('out/significant_interaction_moi_boxplots/')

# read in significant interaction models
significant_results <- get_significant_results()

# read enhancer-guide mapping table
enhancer_guide <- read_enhancer_guide_table()

# read in guide efficiency information
guide_info <- read_guide_efficiency_info()

# read in guide assignment matrix
guide_matrix <- read_guide_matrix()

# read in cell-level covariates
covariates <- read_covariates()

# iterate through significant interactions and plot boxplots
for (i in 1:nrow(significant_results)) {

  # get name of enhancers and gene
  enhancer_1 <- significant_results[i, 'enhancer.1.list']
  enhancer_2 <- significant_results[i, 'enhancer.2.list']
  gene <- significant_results[i, 'gene.list']

  # get guide-efficiency adjusted perturbation vectors for each enhancer
  enh_1_perturbation <- get_enhancer_perturbation(enhancer_1, enhancer_guide,
                                                  guide_info, guide_matrix)
  enh_2_perturbation <- get_enhancer_perturbation(enhancer_2, enhancer_guide,
                                                  guide_info, guide_matrix)

  # create data frame for modeling
  model_df <- cbind(
    enh_1_perturbation,
    enh_2_perturbation,
    covariates
  )

  # get double perturbation MOI
  enh_1_perturbed <- enh_1_perturbation > 0
  enh_2_perturbed <- enh_2_perturbation > 0
  double_perturb <- enh_1_perturbed & enh_2_perturbed
  
  # plot boxplots
  model_df$double_perturb <- double_perturb

  plot <- ggplot(model_df, aes(x = double_perturb, y = guide_count, fill = double_perturb)) +
    geom_boxplot(outlier.shape = 1) +
    theme_classic() +
    xlab('') +
    scale_x_discrete(
      labels = c('All Other', 'Double Perturbation')
    ) +
    scale_fill_manual(
      breaks = c('FALSE', 'TRUE'),
      values = c('darkgray', 'deepskyblue')
    ) + 
    ylab('Guide Count') +
    theme(
      axis.line = element_line(linewidth = 1),
      axis.title.x = element_text(size = 24, color = 'black'),
      axis.title.y = element_text(size = 24, color = 'black'),
      axis.text.x = element_text(size = 16, color = 'black', family = 'Helvetica'),
      axis.text.y = element_text(size = 20, color = 'black', family = 'Helvetica'),
      axis.ticks = element_line(color = 'black', linewidth = 1),
      axis.ticks.length = unit(2, 'mm'),
      plot.margin = rep(unit(10, 'mm'), 4),
      legend.position = "none"
    )


  ggsave(
    plot = plot,
    filename = paste0(
      'out/significant_interaction_moi_boxplots/',
      enhancer_1,
      '_',
      enhancer_2,
      '_',
      gene,
      '.pdf'
    ),
    device = 'pdf',
    unit = 'in',
    width = 5,
    height = 7
  )

}
