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

# compute N50 per gene and 50% counts
gene_umi_counts_n50 <- gene_umi_counts / 2
n50_per_gene <- rep(NA, nrow(expr_matrix))

for (i in 1:5) {

  # get counts for gene 
  gene_counts <- expr_matrix[i, ]

  # sort in descending order
  gene_counts <- sort(gene_counts, decreasing = TRUE)

  # counter number of cells and number of UMIs
  cell_counter <- 0
  num_umis <- 0

  # compute number of cells required to reach 50% of counts
  for (j in 1:length(gene_counts)) {
    cell_counter <- cell_counter + 1
    num_umis <- num_umis + gene_counts[j]

    if (num_umis >= gene_umi_counts_n50[i]) {
      break
    }
  }

  # save to output vector
  n50_per_gene[i] <- cell_counter

}

# create output dataframe
output_df <- data.frame(cbind(prop_umis_metric, n50_per_gene))

# save metric to output file for plotting/visualization
write.csv(
  output_df,
  'data/experimental/processed/umi_cells_n50_metrics.csv',
  quote = FALSE
)
