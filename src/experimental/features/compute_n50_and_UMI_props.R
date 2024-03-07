# There is an effect in the data where a few cells have really high expression
# of a gene. This becomes particularly relevant for the enhancer-enhancer
# interaction models, since a few cells in the double perturbation vector have
# large counts, which leads to large positive interaction coefficients. This
# script generates a UMI proportion based metric for quantifying this effect.
#
# Author: Karthik Guruvayurappan

library(rhdf5)
source('src/experimental/models/analysis_helpers.R')

# read in gene expression matrix
h5_name <- 'data/experimental/processed/gasperini_data.h5'
expr_matrix <- read_expr_matrix(h5_name)

# compute total UMIs per gene
gene_umi_counts <- rowSums(expr_matrix)

# compute number of cells with a UMI per gene
cells_with_umi_counts <- rowSums(expr_matrix > 0)

# take log10
prop_umis_metric <- log10(gene_umi_counts / cells_with_umi_counts)

# save metric to output file for plotting/visualization
write.csv(
  prop_umis_metric,
  'data/experimental/processed/umi_cells_metric.csv',
  quote = FALSE
)
