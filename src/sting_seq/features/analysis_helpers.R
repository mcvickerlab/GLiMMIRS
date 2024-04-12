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


# function to read in and preprocess STING-seq guide matrix
read_guide_matrix <- function() {

  # read in guide matrix using Seurat
  guide_matrix <- Seurat::Read10X(
    data.dir = 'data/sting_seq/raw/guide_matrix'
  )
  guide_matrix <- guide_matrix[[2]]

  # remove '-1' from end of each cell barcode
  colnames(guide_matrix) <- substr(
    colnames(guide_matrix),
    1,
    nchar(colnames(guide_matrix)) - 2
  )

  # read in guide RNA UMI thresholds
  umi_thresholds <- read.csv(
    'data/sting_seq/raw/guide_matrix/umi_thresholds.csv.gz'
  )
  rownames(umi_thresholds) <- umi_thresholds$Protospacer

  # threshold guide RNA assignments by UMI threshold
  for (i in 1:nrow(guide_matrix)) {
    current_guide <- rownames(guide_matrix)[i]
    current_guide_threshold <- umi_thresholds[current_guide, ]$UMI.threshold
    guide_matrix[i, ] <- as.numeric(
      guide_matrix[i, ] >= current_guide_threshold
    )
  }

  # read in supplementary table S3C
  cell_covariates <- read.csv(
    'data/sting_seq/raw/sting_seq_suppl_table_s3c.csv'
  )

  # filter for cells that were included in the STING-seq analysis
  guide_matrix<- guide_matrix[
    , 
    colnames(guide_matrix) %in% cell_covariates$Cell.Barcode
  ]

  return (guide_matrix)
}


# function to read in and process STING-seq covariates
read_covariates <- function() {

  # read in supplementary table S3C
  cell_covariates <- read.csv(
    'data/sting_seq/raw/sting_seq_suppl_table_s3c.csv'
  )

  # filter to only include cells in v1 experiment
  cell_covariates <- cell_covariates[cell_covariates$Library == 'v1', ]

  # get covariates used for GLiMMIRS
  percent_mito <- cell_covariates$Mitochondrial.Percentage
  n_umis <- cell_covariates$Total.Gene.Expression.UMIs
  scaling_factors <- n_umis / 1e6
  grna_counts <- cell_covariates$Unique.gRNAs..after.adaptive.thresholding

  # read in cell cycle scores computed using Seurat
  s_scores <- read.csv('data/sting_seq/interim/cell_cycle_s_scores.csv')
  g2m_scores <- read.csv('data/sting_seq/interim/cell_cycle_g2m_scores.csv')

  s_scores <- s_scores$S.Score
  g2m_scores <- g2m_scores$G2M.Score

  covariates <- data.frame(cbind(
    percent_mito,
    scaling_factors,
    grna_counts,
    s_scores,
    g2m_scores
  ))

  return (covariates)
}


# function to read in SNP-guide mapping table
read_snp_guide <- function() {

  # read in supplementary table S1E (contains PTPRC SNPs)
  ptprc_snps <- read.csv(
    'data/sting_seq/raw/sting_seq_suppl_table_s1e.csv'
  )
  ptprc_snps <- ptprc_snps$rs.ID

  # read in supplementary table S3A (contains guides)
  guide_table <- read.csv(
    'data/sting_seq/raw/sting_seq_suppl_table_s3a.csv'
  )

  # filter for guides of interest
  guide_table <- guide_table[guide_table$gRNA.Library == 'v1', ]
  ptprc_guides <- guide_table[guide_table$Target %in% ptprc_snps, ]

  # filter for necessary columns
  ptprc_guides <- ptprc_guides[, c('Target', 'gRNA.ID')]

  return (ptprc_guides)
}


# function to get SNP perturbation vector from SNP name
get_perturbation_vector <- function(snp_guide_table, snp_name) {

  # get guide names for SNP
  guide_names <- snp_guide_table[snp_guide_table$Target == snp_name, ]$gRNA.ID

  # modify guide name format to match guide matrix
  guide_names <- paste0(
    substr(guide_names, 1, nchar(guide_names) - 2),
    '_',
    substr(guide_names, nchar(guide_names), nchar(guide_names))
  )

  # get rows from guide matrix
  guide_vectors <- data.frame(guide_matrix[guide_names, ])

  # convert to a single SNP perturbation vector
  perturb_vector <- colSums(guide_vectors)
  perturb_vector <- as.numeric(perturb_vector > 0)

  return (perturb_vector)
}

