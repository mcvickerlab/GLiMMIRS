# This script performs bootstrap sampling to generate confidence intervals for
# the significant interactions from the at-scale GLiMMIRS analysis.
#
# Author: Karthik Guruvayurappan

library(stats)
library(rhdf5)
library(MASS)

# useful analysis functions
source('src/gasperini/models/analysis_helpers.R')

# create directory to hold outputs
dir.create('data/gasperini/processed/gasperini_bootstraps/')

# read in significant interactions and filter columns
significant_results <- get_significant_results()
significant_results <- significant_results[
  ,
  c('enhancer.1.list', 'enhancer.2.list', 'gene.list')
]
colnames(significant_results) <- c('enhancer_1', 'enhancer_2', 'gene')

# read in enhancer-guide table
enhancer_guide <- read_enhancer_guide_table()

# read in guide efficiency information
guide_info <- read_guide_efficiency_info()

# read in guide matrix
guide_matrix <- read_guide_matrix()

# read in gene expression matrix
expr_matrix <- read_expr_matrix()

# read in cell-level covariates
covariates <- read_covariates()

for (i in 1:nrow(significant_results)) {

  # define enhancers and gene
  enhancer_1 <- significant_results[i, 'enhancer_1']
  enhancer_2 <- significant_results[i, 'enhancer_2']
  gene <- significant_results[i, 'gene']

  # get perturbation vectors corresponding to enhancers
  enh_1_perturbation <- get_enhancer_perturbation(enhancer_1, enhancer_guide,
                                                  guide_info, guide_matrix)
  enh_2_perturbation <- get_enhancer_perturbation(enhancer_2, enhancer_guide,
                                                  guide_info, guide_matrix)

  # get gene expression counts
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

  # define vectors to hold bootstrap results
  num_bootstraps <- 100
  bootstrap_intercepts <- rep(NA, num_bootstraps)
  bootstrap_enh1_betas <- rep(NA, num_bootstraps)
  bootstrap_enh2_betas <- rep(NA, num_bootstraps)
  bootstrap_interaction_betas <- rep(NA, num_bootstraps)

  for (j in 1:num_bootstraps) {

    # bootstrap cells
    bootstrap_df <- model_df[
      sample(1:nrow(model_df), size = nrow(model_df), replace = TRUE),

    ]

    # fit model
    bootstrap_model <- glm.nb(
      formula = model_formula,
      data = bootstrap_df
    )

    # get coefficients
    model_coeffs <- summary(bootstrap_model)$coefficients
    intercept <- model_coeffs['(Intercept)', 'Estimate']
    enh1_beta <- model_coeffs['enh_1_perturbation', 'Estimate']
    enh2_beta <- model_coeffs['enh_2_perturbation', 'Estimate']
    interaction_beta <- model_coeffs[
      'enh_1_perturbation:enh_2_perturbation',
      'Estimate'
    ]

    # add to output vectors
    bootstrap_intercepts[j] <- intercept
    bootstrap_enh1_betas[j] <- enh1_beta
    bootstrap_enh2_betas[j] <- enh2_beta
    bootstrap_interaction_betas[j] <- interaction_beta
    
  }

  # create output data frame
  output_df <- data.frame(cbind(
    bootstrap_intercepts, 
    bootstrap_enh1_betas, 
    bootstrap_enh2_betas, 
    bootstrap_interaction_betas
  ))

  # write to output file
  write.csv(
    output_df,
    paste0(
      'data/gasperini/processed/gasperini_bootstraps/',
      enhancer_1,
      '_',
      enhancer_2,
      '_',
      gene,
      '_bootstrap_results.csv'
    ),
    row.names = FALSE
  )
}
