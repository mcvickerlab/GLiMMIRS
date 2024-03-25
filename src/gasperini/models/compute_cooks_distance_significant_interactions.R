# This script computes the Cook's distances for all cells in the
# enhancer-enhancer-gene sets that had significant interaction terms from the
# genome-wide analysis.
#
# Author: Karthik Guruvayurappan

library(rhdf5)
library(MASS)
library(stats)

# source helper functions for analysis
source('src/gasperini/models/analysis_helpers.R')

# create output directory
dir.create('data/gasperini/processed/significant_interaction_cooks_distances/')

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

# read in gene expression matrix
expr_matrix <- read_expr_matrix()

# iterate through significant interactions and re-fit GLiMMIRS models
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

  # get gene counts
  gene_counts <- expr_matrix[gene, ]

  # create data frame for modeling
  model_df <- cbind(
    enh_1_perturbation,
    enh_2_perturbation,
    gene_counts,
    covariates
  )

  # fit model
  model_formula <- as.formula(paste0(
    'gene_counts ~ ',
    'enh_1_perturbation * ',
    'enh_2_perturbation + ',
    'prep_batch + ',
    'guide_count + ',
    'percent.mito + ',
    's.score + ',
    'g2m.score + ',
    'offset(log(scaling.factor))'
  ))

  model <- glm.nb(
    formula = model_formula,
    data = model_df
  )

  # compute Cook's distance for each cell
  cooks_distances <- cooks.distance(model)

  # merge data frames together
  model_df <- model_df[complete.cases(model_df), ]
  cooks_df <- data.frame(cbind(model_df, cooks_distances))

  # add perturbation type column
  cooks_df$perturbation_type <- 'No Perturbation'

  enh_1_perturbed <- cooks_df$enh_1_perturbation > 0
  cooks_df$perturbation_type[enh_1_perturbed] <- 'Enhancer 1'

  enh_2_perturbed <- cooks_df$enh_2_perturbation > 0
  cooks_df$perturbation_type[enh_2_perturbed] <- 'Enhancer 2'

  both_perturbed <- (enh_1_perturbed & enh_2_perturbed)
  cooks_df$perturbation_type[both_perturbed] <- 'E1 + E2'

  # subset to necessary columns
  cooks_df <- cooks_df[, c('perturbation_type', 'cooks_distances')]

  # write to output CSV file
  write.csv(
    cooks_df,
    paste0(
      'data/gasperini/processed/significant_interaction_cooks_distances/',
      enhancer_1,
      '_',
      enhancer_2,
      '_',
      gene,
      '.csv'
    ),
    row.names = FALSE
  )

}
