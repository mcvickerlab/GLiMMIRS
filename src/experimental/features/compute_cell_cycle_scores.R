-# This script computes the cell cycle scores for the at-scale screen using the
# Seurat single-cell RNA-sequencing analysis package.
#
# Author: Karthik Guruvayurappan

library(Matrix)
library(biomaRt)
library(Seurat)
library(ggplot2)

# load in UMI count (expression) matrix from at-scale screen
data.path <- 'data/experimental/raw/'
expression.matrix <- readMM(
    paste0(data.path, 'GSE120861_at_scale_screen.exprs.mtx.gz')
)

# convert expression matrix from matrix to data frame format
expression.matrix <- as.data.frame(expression.matrix)

# read in column names and add to expression matrix
cell.barcodes <- read.delim(
    paste0(data.path, 'GSE120861_at_scale_screen.cells.txt.gz'),
    header = FALSE
)
cell.barcodes <- cell.barcodes$V1
colnames(expression.matrix) <- cell.barcodes

# read in row names and add to expression matrix
genes <- read.delim(
    paste0(data.path, 'GSE120861_at_scale_screen.genes.txt.gz'),
    header = FALSE
)
genes <- genes$V1
rownames(expression.matrix) <- genes

# use BioMart to convert Ensembl gene IDs to HGNC symbols
mart <- useMart('ensembl')
mart <- useDataset('hsapiens_gene_ensembl', mart)
gene.symbols <- getBM(
    filters = 'ensembl_gene_id',
    attributes = c('ensembl_gene_id', 'hgnc_symbol'),
    values = genes,
    mart = mart,
)

# merge with original gene symbols and retain gene ordering
gene.symbols <- merge(
    data.frame(genes),
    gene.symbols,
    all.x = TRUE,
    by.x = 'genes',
    by.y = 'ensembl_gene_id',
    sort = FALSE
)

# merge again to avoid NAs at bottom of data frame
gene.symbols <- merge(
    data.frame(genes),
    gene.symbols,
    by = 'genes',
    sort = FALSE
)

# impute Ensembl Gene ID for genes without HGNC symbol
gene.symbols[is.na(gene.symbols)] <- ''

for (i in 1:nrow(gene.symbols)) {
    if (gene.symbols[i, 'hgnc_symbol'] == '') {
        gene.symbols[i, 'hgnc_symbol'] <- gene.symbols[i, 'genes']
    }
}

# set gene names of expression matrix to HGNC symbols
rownames(expression.matrix) <- gene.symbols[, 'hgnc_symbol']

# use predefined S and G2M gene sets (from Seurat)
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

# run Seurat initialization steps for expression matrix
gene.expression <- CreateSeuratObject(
    counts = expression.matrix
)
gene.expression <- NormalizeData(gene.expression)
gene.expression <- FindVariableFeatures(
    gene.expression,
    selection.method  = "vst"
)
gene.expression <- ScaleData(
    gene.expression,
    features = rownames(gene.expression)
)

# run Seurat cell cycle scoring method
gene.expression <- CellCycleScoring(
    gene.expression,
    s.features = s.genes,
    g2m.features = g2m.genes,
    set.ident = TRUE
)

# run PCA to show cell separation based on cell cycle phase
gene.expression <- RunPCA(gene.expression, features = c(s.genes, g2m.genes))
cell.cycle.pca <- DimPlot(gene.expression, raster = FALSE)

# save plot
ggsave(
    filename = 'out/cell_cycle_pca.png',
    plot = cell.cycle.pca,
    device = 'png',
    height = 7,
    width = 7,
    units = 'in'
)

# write scores to CSV files
s.scores <- gene.expression[[]]['S.Score']
g2m.scores <- gene.expression[[]]['G2M.Score']

write.csv(s.scores, 'data/experimental/interim/cell_cycle_s_scores.csv')
write.csv(g2m.scores, 'data/experimental/interim/cell_cycle_g2m_scores.csv')
