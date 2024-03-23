# To determine p-values associated with the significant interaction
# coefficients that were sensitive to outliers, we ran 100 permutation tests
# for each enhancer pair and recorded the interaction coefficients. This
# code plots the distribution of those coefficients, along with where the
# observed interaction coefficient falls on that distribution.
#
# Author: Karthik Guruvayurappan

source('src/gasperini/models/analysis_helpers.R')
source('src/gasperini/visualization/plotting_helpers.R')

library(ggplot2)

# create directory to hold output files
dir.create('out/adaptive_permutation_histograms/')

# get significant interactions from at-scale analysis
significant_interactions <- get_significant_results()

# iterate through significant interactions and plot permutation histograms
for (i in 1:nrow(significant_interactions)) {

  # get enhancer 1, enhancer 2, and gene name
  enhancer_1 <- significant_interactions[i, 'enhancer.1.list']
  enhancer_2 <- significant_interactions[i, 'enhancer.2.list']
  gene <- significant_interactions[i, 'gene.list']
  pair_name <- concat_enhancer_pair_name(enhancer_1, enhancer_2, gene)

  # store value of observed interaction coefficient
  observed_interaction <- significant_interactions[i, 'interaction.effects']

  # store filename for permutations and read in permutation results
  permutation_results_filename <- paste0(
    'data/experimental/processed/adaptive_permutation_test_results/',
    pair_name,
    '_adaptive_permutations.csv'
  )

  # read in file containing permutation test results
  permutation_results <- read.csv(permutation_results_filename)

  # plot histograms of permutation coefficients
  plot <- plot_histogram(permutation_results, interaction.effects) +
    geom_vline(
      xintercept = observed_interaction,
      color = 'red',
      linetype = 'dashed'
    ) +
    xlab('Interaction Coefficient') +
    ylab('Count') +
    ggtitle(paste(enhancer_1, enhancer_2, gene)) +
    theme(
      axis.title.x = element_text(size = 14, color = 'black'),
      axis.title.y = element_text(size = 14, color = 'black'),
      axis.text = element_text(size = 10, color = 'black'),
      plot.title = element_text(size = 8.5, color = 'black', hjust = 0.5)
    )

  # save to output files
  ggsave(
    filename = paste0(
      'out/adaptive_permutation_histograms/',
      pair_name,
      '.png'
    ),
    device = 'png',
    plot = plot
  )
}
