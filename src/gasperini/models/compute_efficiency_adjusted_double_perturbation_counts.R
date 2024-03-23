# This script re-runs the GLiMMIRS model on the significant interactions and
# records the Cook's distance associated with each cell, along with its
# perturbation.
#
# Author: Karthik Guruvayurappan

source('src/gasperini/models/analysis_helpers.R')

library(rhdf5)
library(MASS)

# read in enhancer pairs
enhancer_pairs <- read.csv(
    'data/gasperini/processed/perturbation_counts.csv'
)

# read in gene names (to filter pairs)
h5_name <- 'data/gasperini/processed/gasperini_data.h5'
genes <- h5read(
    h5_name,
    'expr/gene_names'
)
enhancer_pairs <- enhancer_pairs[enhancer_pairs$gene %in% genes, ]

# filter to >= 10 cell threshold
enhancer_pairs <- enhancer_pairs[enhancer_pairs$double.count.list >= 10, ]
colnames(enhancer_pairs) <- c('enhancer_1', 'enhancer_2', 'gene')
enhancer_pairs <- enhancer_pairs[, c('enhancer_1', 'enhancer_2', 'gene')]

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

# initialize empty vectors to hold "true" double perturbation counts
enhancer_1_list <- rep(NA, nrow(enhancer_pairs))
enhancer_2_list <- rep(NA, nrow(enhancer_pairs))
gene_list <- rep(NA, nrow(enhancer_pairs))
double_perturbation_counts <- rep(NA, nrow(enhancer_pairs))

for (i in 1:10) {

  # get name of enhancers and gene
  enhancer_1 <- enhancer_pairs[i, 'enhancer_1']
  enhancer_2 <- enhancer_pairs[i, 'enhancer_2']
  gene <- enhancer_pairs[i, 'gene']

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

  # record double perturbation count info
  double_perturbation_count <- sum(
    model_df$enh_1_perturbation * model_df$enh_2_perturbation > 0
  )
  enhancer_1_list[i] <- enhancer_1
  enhancer_2_list[i] <- enhancer_2
  gene_list[i] <- gene
  double_perturbation_counts[i] <- double_perturbation_count                        
}

# write double perturbation counts to output file
output_df <- data.frame(cbind(
  enhancer_1_list,
  enhancer_2_list,
  gene_list,
  double_perturbation_counts
))

write.csv(
  output_df, 
  paste0(
    'data/gasperini/processed/',
    'enhancer_pair_efficiency_adjusted_double_perturb_counts.csv'
  ),
  row.names = FALSE
)
