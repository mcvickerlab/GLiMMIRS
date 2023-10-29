# This script computes cell cycle scores for the STING-seq data.
#
# Author: Karthik Guruvayurappan

library(Matrix) # for Matrix reading
library(Seurat) # for single-cell analysis
library(ggplot2) # for plotting

# read in STING-seq expression matrix (for all lanes)
rna.a <- Read10X(
    '/iblm/netapp/data1/external/Morris_2023_STING_seq/cDNA/A'
)
rna.b <- Read10X(
    '/iblm/netapp/data1/external/Morris_2023_STING_seq/cDNA/B'
)
rna.c <- Read10X(
    '/iblm/netapp/data1/external/Morris_2023_STING_seq/cDNA/C'
)
rna.d <- Read10X(
    '/iblm/netapp/data1/external/Morris_2023_STING_seq/cDNA/D'
)

# merge all of the matrices together column-wise
rna <- cbind(rna.a, rna.b, rna.c, rna.d)

# read in supplementary table S3C
table.s3c <- read.csv(
    '/iblm/netapp/data1/external/Morris_2023_STING_seq/sting_seq_suppl_table_s3c.csv'
)


# use predefined S and G2M gene sets (from Seurat)
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

# run Seurat initialization steps for expression matrix
rna <- CreateSeuratObject(
    counts = rna
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
    s.features = s.genes,
    g2m.features = g2m.genes,
    set.ident = TRUE
)

# run PCA to show cell separation based on cell cycle phase
rna <- RunPCA(rna, features = c(s.genes, g2m.genes))
cell.cycle.pca <- DimPlot(rna, raster = FALSE)

# save plot to output file
ggsave(
    filename = '/iblm/netapp/home/karthik/GLiMMIRS/out/23_10_29_cell_cycle_pca_laneA.png',
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





