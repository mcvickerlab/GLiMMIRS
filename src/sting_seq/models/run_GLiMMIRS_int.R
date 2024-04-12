# This script runs GLiMMIRS-int on all pairwise combinations of SNPs for the
# PTPRC locus to see if there are interactions between regulatory sequences.
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

# get unique SNPs for PTPRC
ptprc_snps <- unique(snp_guide_table$Target)

# generate pairwise combinations of SNPs
ptprc_pairs <- data.frame(t(combn(ptprc_snps, 2)))
colnames(ptprc_pairs) <- c('snp_1', 'snp_2')

# create vectors to hold output results
snp_1_estimates <- rep(NA, nrow(ptprc_pairs))
snp_1_pvalues <- rep(NA, nrow(ptprc_pairs))

snp_2_estimates <- rep(NA, nrow(ptprc_pairs))
snp_2_pvalues <- rep(NA, nrow(ptprc_pairs))

interaction_estimates <- rep(NA, nrow(ptprc_pairs))
interaction_pvalues <- rep(NA, nrow(ptprc_pairs))

# iterate through SNP pairs and run GLiMMIRS-int
for (i in 1:nrow(ptprc_pairs)) {

  # get SNP names
  snp_1 <- ptprc_pairs$snp_1[i]
  snp_2 <- ptprc_pairs$snp_2[i]

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

  # save estimates to output vectors
  mdl_results <- summary(mdl)$coefficients
  
  snp_1_estimates[i] <- mdl_results['snp_1_perturb', 'Estimate']
  snp_1_pvalues[i] <- mdl_results['snp_1_perturb', 'Pr(>|z|)']

  snp_2_estimates[i] <- mdl_results['snp_2_perturb', 'Estimate']
  snp_2_pvalues[i] <- mdl_results['snp_2_perturb', 'Pr(>|z|)']

  interaction_estimates[i] <- mdl_results[
    'snp_1_perturb:snp_2_perturb',
    'Estimate'
  ]
  interaction_pvalues[i] <- mdl_results[
    'snp_1_perturb:snp_2_perturb',
    'Pr(>|z|)'
  ]
}

# create results data frame
results_df <- data.frame(cbind(
  ptprc_pairs,
  snp_1_estimates,
  snp_1_pvalues,
  snp_2_estimates,
  snp_2_pvalues,
  interaction_estimates,
  interaction_pvalues
))

# perform FDR and Bonferroni correction on interaction pvalues
results_df$fdr_int_pvalues <- p.adjust(interaction_pvalues, method = 'fdr')
results_df$bonf_int_pvalues <- p.adjust(
  interaction_pvalues, 
  method = 'bonferroni'
)

# write output file to CSV
write.csv(
  results_df,
  'data/sting_seq/processed/GLiMMIRS_int_ptprc.csv',
  row.names = FALSE
)
