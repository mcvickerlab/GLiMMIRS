# This script plots a histogram of the Cook's distances for the double
# perturbation cells where a significant interaction is observed in the
# GLiMMIRS model.
#
# Author: Karthik Guruvayurappan

library(ggplot2)

# source helper analysis functions
source('src/gasperini/models/analysis_helpers.R')

# create directory to hold output plots
dir.create('out/significant_interaction_cooks_distances/')

# get significant interactions
significant_results <- get_significant_results()

# record summary statistics for cook's distances
enhancer_1_list <- rep(NA, nrow(significant_results))
enhancer_2_list <- rep(NA, nrow(significant_results))
gene_list <- rep(NA, nrow(significant_results))
mean_cooks_distances <- rep(NA, nrow(significant_results))
max_cooks_distances <- rep(NA, nrow(significant_results))
mean_double_perturb_cooks_distances <- rep(NA, nrow(significant_results))
discard_cooks_distance <- rep(NA, nrow(significant_results))

# iterate through significant interactions
for (i in 1:nrow(significant_results)) {

  # read in perturbation-Cook's distance data frame
  enhancer_1 <- significant_results[i, 'enhancer.1.list']
  enhancer_2 <- significant_results[i, 'enhancer.2.list']
  gene <- significant_results[i, 'gene.list']
  
  cooks_filename <- paste0(
    'data/gasperini/processed/significant_interaction_cooks_distances/',
    enhancer_1, 
    '_',
    enhancer_2,
    '_',
    gene,
    '.csv'
  )

  cooks_distances <- read.csv(cooks_filename)

  # print top 50 cells for Cook's distance
  cooks_distances <- cooks_distances[
    order(cooks_distances$cooks_distances, decreasing = TRUE),
    
  ]

  # compute summary statistics and add to output
  mean_cooks_distance <- mean(cooks_distances$cooks_distances)
  max_cooks_distance <- max(cooks_distances$cooks_distances)
  mean_cooks_distances[i] <- mean_cooks_distance
  max_cooks_distances[i] <- max_cooks_distance


  # plot histogram for cells with the double perturbation cells
  double_perturb_cells <- cooks_distances[
    cooks_distances$perturbation_type == 'E1 + E2',
    
  ]

  # compute mean of double perturb cell cook's distances
  mean_double_perturb_cooks_distances[i] <- mean(
    double_perturb_cells$cooks_distances
  )

  # compute whether enhancer pair would be discarded or not
  discard_pair <- (
    (max_cooks_distance / mean(double_perturb_cells$cooks_distances)) > 5
  )
  discard_cooks_distance[i] <- discard_pair

  # plot double perturb cell cook's distances
  plot <- ggplot(double_perturb_cells, aes(x = cooks_distances)) +
    geom_histogram(fill = 'darkgray', color = 'darkgray') +
    theme_classic() +
    scale_x_continuous(expand = c(0.02, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    xlab("Cook's Distance") + 
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
    filename = paste0(
      'out/significant_interaction_cooks_distances/',
      enhancer_1,
      '_',
      enhancer_2,
      '_',
      gene,
      '.pdf'
    ),
    plot = plot,
    device = 'pdf'
  )

  # save enhancer names and gene name to output list
  enhancer_1_list[i] <- enhancer_1
  enhancer_2_list[i] <- enhancer_2
  gene_list[i] <- gene
}

# create output dataframe of summary statistics
summary_df <- data.frame(cbind(
  enhancer_1_list,
  enhancer_2_list,
  gene_list,
  mean_cooks_distances,
  max_cooks_distances,
  mean_double_perturb_cooks_distances,
  discard_cooks_distance
))

# write summary statistics to output file
write.csv(
  summary_df,
  'data/gasperini/processed/interaction_cooks_distance_summary_stats.csv',
  row.names = FALSE
)
