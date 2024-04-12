# This script computes the Cook's distances associated with each cell for the
# enhancer-enhancer-gene sets with significant interaction terms.
#
# Author: Karthik Guruvayurappan

library(Matrix)
library(Seurat)
library(MASS)
library(dplyr)

# read in useful helper functions for STING-seq analysis
source('src/sting_seq/features/analysis_helpers.R')

# read in expression matrix
expr_matrix <- read_expr_matrix()

# read in guide matrix
guide_matrix <- read_guide_matrix()

# read in cell-level covariates
covariates <- read_covariates()

# read in SNP-guide table
snp_guide_table <- read_snp_guide()

# read in significant interaction results
results <- read.csv('data/sting_seq/processed/GLiMMIRS_int_ptprc.csv')
significant_results <- results[results$fdr_int_pvalues < 0.1, ]

# iterate through significant results
for (i in 1:nrow(significant_results)) {

  # get SNP names
  snp_1 <- significant_results[i, 'snp_1']
  snp_2 <- significant_results[i, 'snp_2']

  # get perturbation vectors for each SNP
  snp_1_perturb <- get_perturbation_vector(
    snp_guide_table,
    snp_1
  )

  snp_2_perturb <- get_perturbation_vector(
    snp_guide_table,
    snp_2
  )

  # get PTPRC expression
  ptprc <- expr_matrix['PTPRC', ]

  # create dataframe for modeling
  model_df <- data.frame(cbind(
    snp_1_perturb,
    snp_2_perturb,
    ptprc,
    covariates
  ))

  # fit negative binomial GLM model
  mdl <- glm.nb(
    ptprc ~ snp_1_perturb * snp_2_perturb + percent_mito + grna_counts + s_scores + g2m_scores + offset(log(scaling_factors)), data = model_df
  )

  # get Cook's distances
  cooks_distances <- cooks.distance(mdl)

  # make data frame for Cook's distance
  cooks_df <- data.frame(cbind(model_df, cooks_distances))
  cooks_df$perturbation <- 'No Perturbation'

  enh_1_perturbed <- (snp_1_perturb > 0)
  cooks_df$perturbation[enh_1_perturbed] <- 'Enhancer 1'

  enh_2_perturbed <- (snp_2_perturb > 0)
  cooks_df$perturbation[enh_2_perturbed] <- 'Enhancer 2'

  both_perturbed <- (enh_1_perturbed & enh_2_perturbed)
  cooks_df$perturbation[both_perturbed] <- 'E1 + E2'

  # subset to necessary columns
  cooks_df <- cooks_df[, c('perturbation', 'cooks_distances')]

  # write Cook's distances to output files
  write.csv(
    cooks_df,
    paste0(
      'data/sting_seq/processed/',
      snp_1,
      '_',
      snp_2,
      '_cooks_distances.csv'
    ),
    row.names = FALSE
  )
}



