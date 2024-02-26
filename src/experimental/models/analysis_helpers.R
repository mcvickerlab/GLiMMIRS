# This file contains numerous helper functions for GLiMMIRS analysis.
#
# Author: Karthik Guruvayurappan

#' Read in at-scale analysis and filter for significant results
#'
get_significant_results <- function() {

  # read in results from at-scale enhancer pair analysis
  model.results <- data.frame()

  for (i in 1:32) {
    batch.file <- paste0(
      'data/experimental/processed/enhancer_pairs_at_scale_',
      i,
      '.csv'
    )
    batch.results <- read.csv(batch.file)
    model.results <- rbind(model.results, batch.results)
  }

  # filter for cases with valid guide efficiencies
  model.results <- model.results[stats::complete.cases(model.results), ]

  # compute FDR-adjusted p-values and filter for significant results
  model.results$adj.interaction.pvalues <- p.adjust(
      model.results$interaction.pvalues,
      method = 'fdr'
  )
  significant.results <- model.results[
      model.results$adj.interaction.pvalues < 0.1,

  ]

  significant.results
}

