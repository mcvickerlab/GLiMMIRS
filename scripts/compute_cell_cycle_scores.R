# This script computes the cell cycle scores for the at-scale screen using the
# Seurat single-cell RNA-sequencing analysis package.
#
# Author: Karthik Guruvayurappan

library('Seurat')
library('Matrix')
library('biomaRt')
library(ggplot2)

# load in UMI count (expression) matrix from at-scale screen
expression.matrix <- readMM('/iblm/netapp/data1/external/Gasperini2019/suppl/GSE120861_at_scale_screen.exprs.mtx')

# convert expression matrix from matrix to data frame format
expression.matrix <- as.data.frame(expression.matrix)

# read in column names and add to expression matrix
cell.barcodes <- read.delim('/iblm/netapp/data1/external/Gasperini2019/suppl/GSE120861_at_scale_screen.cells.txt', header = FALSE)
cell.barcodes <- cell.barcodes$V1
colnames(expression.matrix) <- cell.barcodes

# read in row names and add to expression matrix
genes <- read.delim('/iblm/netapp/data1/external/Gasperini2019/suppl/GSE120861_at_scale_screen.genes.txt', header = FALSE)
genes <- genes$V1
rownames(expression.matrix) <- genes

# code snippet adapted from: https://stackoverflow.com/questions/28543517/how-can-i-convert-ensembl-id-to-gene-symbol-in-r
mart <- useMart('ensembl')
mart <- useDataset("hsapiens_gene_ensembl", mart)
gene.symbols <- getBM(
    filters = "ensembl_gene_id",
    attributes = c("ensembl_gene_id", "hgnc_symbol"),
    values = genes,
    mart = mart,
)

# merge with original gene symbols to retain gene ordering and add back genes that had no match
gene.symbols <- merge(data.frame(genes), gene.symbols, all.x = TRUE, by.x = 'genes', by.y = 'ensembl_gene_id', sort = FALSE)

gene.symbols <- merge(data.frame(genes), gene.symbols, by = 'genes', sort = FALSE)

gene.symbols[is.na(gene.symbols)] <- ''

# set gene names to HGNC symbols
for (i in 1:nrow(gene.symbols)) {
    if (gene.symbols[i, 'hgnc_symbol'] == '') {
        gene.symbols[i, 'hgnc_symbol'] <- gene.symbols[i, 'genes']
    }
}
rownames(expression.matrix) <- gene.symbols[, 'hgnc_symbol']

# use predefined S and G2M gene sets
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

# run Seurat initialization steps
gene.expression <- CreateSeuratObject(counts = expression.matrix)
gene.expression <- NormalizeData(gene.expression)
gene.expression <- FindVariableFeatures(gene.expression, selection.method  = "vst")
gene.expression <- ScaleData(gene.expression, features = rownames(gene.expression))
gene.expression <- CellCycleScoring(gene.expression, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

# Running a PCA on cell cycle genes reveals, unsurprisingly, that cells separate entirely by
# phase
gene.expression <- RunPCA(gene.expression, features = c(s.genes, g2m.genes))
cell.cycle.pca <- DimPlot(gene.expression)

# save plot
ggsave(
    filename = '/iblm/netapp/home/karthik/GLiMMIRS/plots/23_03_26_cell_cycle_pca.pdf',
    plot = cell.cycle.pca,
    device = 'pdf'
)

# write scores to CSV files
s.scores <- gene.expression[[]]['S.Score']
g2m.scores <- gene.expression[[]]['G2M.Score']

write.csv(s.scores, '/iblm/netapp/home/karthik/GLiMMIRS/gasperini_data/s_scores.csv')
write.csv(g2m.scores, '/iblm/netapp/home/karthik/GLiMMIRS/gasperini_data/g2m_scores.csv')
