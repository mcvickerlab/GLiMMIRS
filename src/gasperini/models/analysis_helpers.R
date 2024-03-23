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
      'data/gasperini/processed/enhancer_pairs_at_scale_',
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

# define name of h5 file as variable
h5_name <- 'data/gasperini/processed/gasperini_data.h5'

# Read in enhancer-guide mapping table
read_enhancer_guide_table <- function() {
  
  # read in enhancer-guide table
  enhancer_guide <- h5read(
      h5_name,
      'enhancer_guide'
  )

  return (enhancer_guide)
}


# read in guide efficiency info table
read_guide_efficiency_info <- function() {
  guide_info <- h5read(
      h5_name,
      'grna/guide_info'
  )
  return(guide_info)
}


# read in guide matrix
read_guide_matrix <- function() {
  # read in guide matrix
  guide_matrix <- h5read(
      h5_name,
      'grna/guide_matrix'
  )
  guide_names <- h5read(
      h5_name,
      'grna/guide_names'
  )
  rownames(guide_matrix) <- guide_names
  barcodes <- h5read(
      h5_name,
      'grna/cell_barcodes'
  )
  colnames(guide_matrix) <- barcodes

  return(guide_matrix)
}


# read in cell-level covariates
read_covariates <- function() {
  covariates <- h5read(
      h5_name,
      'expr/cell_covariates'
  )
  covariates$scaling.factor <- as.vector(covariates$scaling.factor)

  return(covariates)
}


# read in expression matrix
read_expr_matrix <- function() {
  # read in counts matrix
  expr_matrix <- h5read(
      h5_name,
      'expr/expr_matrix'
  )
  genes <- h5read(
      h5_name,
      'expr/gene_names'
  )
  rownames(expr_matrix) <- genes
  barcodes <- h5read(
      h5_name,
      'expr/cell_barcodes'
  )
  colnames(expr_matrix) <- barcodes

  # add pseudocount to count data
  pseudocount <- 0.01
  expr_matrix <- expr_matrix + pseudocount

  return(expr_matrix)
}


# get guide efficiency-adjusted enhancer perturbation vector
get_enhancer_perturbation <- function(enhancer, enhancer_guide, guide_info,
                                      guide_matrix) {
  
  # get spacer sequences for enhancer
  enhancer_spacers <- enhancer_guide[
    enhancer_guide$target.site == enhancer,
    'spacer'
  ]

  # get guide efficiency info for spacers
  enh_efficiencies <- guide_info[guide_info$spacer %in% enhancer_spacers, ]
  enh_efficiencies <- enh_efficiencies[
    ,
    c('spacer', 'Cutting.Efficiency')
  ]
  enh_efficiencies[is.na(enh_efficiencies)] <- 0

  # compute perturbation vector using guide efficiency info
  enh_no_perturbation <- rep(1, ncol(guide_matrix))
  for (i in 1:nrow(enh_efficiencies)) {
      spacer <- enh_efficiencies[i, 'spacer']
      efficiency <- enh_efficiencies[i, 'Cutting.Efficiency']
      spacer_vector <- guide_matrix[spacer, ]
      spacer_perturbation <- spacer_vector * efficiency
      spacer_no_perturbation <- 1 - spacer_perturbation
      enh_no_perturbation <- enh_no_perturbation * spacer_no_perturbation
  }
  enh_perturbation <- 1 - enh_no_perturbation

  return (enh_perturbation)
}


