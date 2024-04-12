# This file contains some helper functions for common code chunks that appear
# in the analysis of the STING-seq data.
#
# Author: Karthik Guruvayurappan

# function to read in STING-seq gene expression matrix
read_expr_matrix <- function() {
  
  # read in expression matrix using Seurat
  expr_matrix <- Seurat::Read10X(
    data.dir = 'data/sting_seq/raw/gene_expression/'
  )

  # remove '-1' string at the end of each cell barcode
  colnames(expr_matrix) <- substr(
    colnames(expr_matrix), 
    1, 
    nchar(colnames(expr_matrix)) - 2
  )

  # read in supplementary table S3C (cell-level covariates)
  cell_covariates <- read.csv(
    'data/sting_seq/raw/sting_seq_suppl_table_s3c.csv'
  )

  # filter for cells that were included in STING-seq analysis
  expr_matrix <- expr_matrix[
    , 
    colnames(expr_matrix) %in% cell_covariates$Cell.Barcode
  ]

  return (expr_matrix)
}

