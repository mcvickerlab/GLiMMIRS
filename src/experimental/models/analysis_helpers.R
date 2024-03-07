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

#' Read in entire at-scale expression matrix
#'
read_expr_matrix <- function(h5_name) {

  # read in expression matrix values
  expr_matrix <- rhdf5::h5read(
    h5_name,
    'expr/expr_matrix'
  )

  # add gene names as rows
  genes <- rhdf5::h5read(
    h5_name,
    'expr/gene_names'
  )
  rownames(expr_matrix) <- genes

  # add cell barcodes as columns
  barcodes <- rhdf5::h5read(
    h5_name,
    'expr/cell_barcodes'
  )
  colnames(expr_matrix) <- barcodes
  
  return(expr_matrix)
}

