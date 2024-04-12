# This script computes cell cycle scores for the STING-seq data.
#
# Author: Karthik Guruvayurappan

library(Matrix)
library(Seurat)
library(ggplot2)

# source analysis helper functions
source('src/sting_seq/features/analysis_helpers.R')

# read in STING-seq expression matrix
expr_matrix <- read_expr_matrix()

# use predefined S and G2M gene sets (from Seurat)
s_genes <- cc.genes$s.genes
g2m_genes <- cc.genes$g2m.genes

# run Seurat initialization steps for expression matrix
rna <- CreateSeuratObject(
  counts = expr_matrix
)
rna <- NormalizeData(rna)
rna <- FindVariableFeatures(
  rna,
  selection.method  = "vst"
)
rna <- ScaleData(
  rna,
  features = rownames(rna)
)

# run Seurat cell cycle scoring method
rna <- CellCycleScoring(
  rna,
  s.features = s_genes,
  g2m.features = g2m_genes,
  set.ident = TRUE
)

# write scores to CSV files
s_scores <- rna[[]]['S.Score']
g2m_scores <- rna[[]]['G2M.Score']

write.csv(
  s_scores,
  'data/sting_seq/interim/cell_cycle_s_scores.csv'
)
write.csv(
  g2m_scores,
  'data/sting_seq/interim/cell_cycle_g2m_scores.csv'
)

