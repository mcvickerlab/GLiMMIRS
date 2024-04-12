# This script runs GLiMMIRS-base to determine enhancer-gene pairs present at
# the PTPRC locus.
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

# create output vectors to hold summary statistics
snps <- rep(NA, length(ptprc_snps))
estimates <- rep(NA, length(ptprc_snps))
pvalues <- rep(NA, length(ptprc_snps))

# run GLiMMIRS-base for each SNP
for (i in 1:length(ptprc_snps)) {

  # get current SNP name
  snp <- ptprc_snps[i]
  snps[i] <- snp

  # get perturbation vector
  snp_perturb <- get_perturbation_vector(
    snp_guide_table,
    snp
  )

  # get PTPRC expression
  ptprc <- expr_matrix['PTPRC', ]

  # create data frame for modeling
  model_df <- data.frame(cbind(
    snp_perturb,
    ptprc,
    covariates
  ))

  # fit negative binomial GLM
  mdl <- glm.nb(
    ptprc ~ snp_perturb + percent_mito + grna_counts + s_scores + g2m_scores + offset(log(scaling_factors)),
    data = model_df
  )

  # get model estimate and pvalue
  estimates[i] <- summary(mdl)$coefficients['snp_perturb', 'Estimate']
  pvalues[i] <- summary(mdl)$coefficients['snp_perturb', 'Pr(>|z|)']
}

# create output data frame 
results_df <- data.frame(cbind(ptprc_snps, estimates, pvalues))

# perform FDR and Bonferroni correction
results_df$fdr_pvalues <- p.adjust(pvalues, method = 'fdr')
results_df$bonferroni_pvalues <- p.adjust(pvalues, method = 'bonferroni')

write.csv(
  results_df,
  'data/sting_seq/processed/GLiMMIRS_base_ptprc_results.csv',
  row.names = FALSE
)


