# This script creates a single HDF5 data structure that contains all of the
# Gasperini data tables and matrices for easy downstream modeling.
#
# Author: Karthik Guruvayurappan

library(Matrix)
library(rhdf5)
library(readxl)

# create an empty h5 file
h5.name <- 'data/experimental/processed/gasperini_data.h5'
h5createFile(h5.name)

# create group for guide and expression matrix
h5createGroup(h5.name, 'grna')
h5createGroup(h5.name, 'expr')

# read in enhancer-gene pairs (for baseline models)
enhancer.gene <- read_excel(
    'data/experimental/raw/suppl_table_2.xlsx',
    sheet = 3
)
enhancer.gene <- data.frame(enhancer.gene)

# data frame should have a shape of 664 x 11

# write to h5 structure
h5write(
    enhancer.gene,
    h5.name,
    'enhancer_gene'
)

# read in smaller subset of enhancer-enhancer pairs
enhancer.enhancer.330 <- read.csv(
    'data/experimental/processed/enhancer_pairs_suppl_table_2.csv'
)
enhancer.enhancer.330$enhancer_1 <- paste0(
    tolower(substring(enhancer.enhancer.330$enhancer_1, 1, 3)),
    substring(
        enhancer.enhancer.330$enhancer_1,
        4,
        nchar(enhancer.enhancer.330$enhancer_1)
    )
)
enhancer.enhancer.330$enhancer_2 <- paste0(
    tolower(substring(enhancer.enhancer.330$enhancer_2, 1, 3)),
    substring(
        enhancer.enhancer.330$enhancer_2,
        4,
        nchar(enhancer.enhancer.330$enhancer_2)
    )
)
enhancer.enhancer.330$enhancer_1 <- sapply(
    enhancer.enhancer.330$enhancer_1,
    FUN = function(x) {
        if (startsWith(x, 'chr')) {
            return (strsplit(x, '_')[[1]][1])
        }
        else {
            return (x)
        }
    }
)
enhancer.enhancer.330$enhancer_2 <- sapply(
    enhancer.enhancer.330$enhancer_2,
    FUN = function(x) {
        if (startsWith(x, 'chr')) {
            return (strsplit(x, '_')[[1]][1])
        }
        else {
            return (x)
        }
    }
)

# write to h5 structure
h5write(
    enhancer.enhancer.330,
    h5.name,
    'enhancer_enhancer_330'
)

# read in at-scale set of enhancer-enhancer pairs
enhancer.enhancer.at.scale <- read.csv(
    'data/experimental/processed/enhancer_pairs_at_scale.csv'
)

# write to h5 structure
h5write(
    enhancer.enhancer.at.scale,
    h5.name,
    'enhancer_enhancer_at_scale'
)

# read in enhancer-guide table
enhancer.guide <- read.table(
    'data/experimental/raw/GSE120861_grna_groups.at_scale.txt.gz'
)
colnames(enhancer.guide) <- c('target.site', 'spacer')
enhancer.guide$target.site <- sapply(
    enhancer.guide$target.site,
    FUN = function(x) {
        if (startsWith(x, 'chr')) {
            return (strsplit(x, '_')[[1]][1])
        }
        else {
            return (x)
        }
    }
)

# write to h5 structure
h5write(
    enhancer.guide,
    h5.name,
    'enhancer_guide'
)

# read in guide-level metadata from GuideScan
guide.info <- read.csv(
    'data/experimental/interim/guidescan_results.csv'
)
guide.info$spacer <- substring(
    guide.info$gRNA,
    1,
    nchar(guide.info$gRNA) - 3
)

# write to h5 structure
h5write(
    guide.info,
    h5.name,
    'grna/guide_info'
)

# read in expression matrix
expr.matrix <- readMM(
    'data/experimental/raw/GSE120861_at_scale_screen.exprs.mtx.gz'
)

# convert to dense matrix representation (for h5)
expr.matrix <- as.matrix(expr.matrix)

# define data dimensionality and chunk size
h5createDataset(
    h5.name,
    'expr/expr_matrix',
    c(nrow(expr.matrix), ncol(expr.matrix)),
    storage.mode = 'integer',
    chunk = c(2000, ncol(expr.matrix))
)

# write to h5 file
h5write(
    expr.matrix,
    h5.name,
    'expr/expr_matrix'
)


# write gene names to h5 structure
genes <- read.table(
    'data/experimental/raw/GSE120861_at_scale_screen.genes.txt.gz'
)
genes <- genes$V1

h5write(
    genes,
    h5.name,
    'expr/gene_names'
)

# write cell barcodes to h5 structure
barcodes <- read.table(
    'data/experimental/raw/GSE120861_at_scale_screen.cells.txt.gz'
)
barcodes <- barcodes$V1

h5write(
    barcodes,
    h5.name,
    'expr/cell_barcodes'
)

# read in phenodata file (for covariates)
covariates <- read.table(
    'data/experimental/raw/GSE120861_at_scale_screen.phenoData.txt.gz'
)
colnames(covariates) <- c(
    'sample',
    'cell',
    'total_umis',
    'size_factor',
    'gene',
    'all_gene',
    'barcode',
    'read_count',
    'umi_count',
    'proportion',
    'guide_count',
    'sample_directory',
    'ko_barcode_file',
    'id',
    'prep_batch',
    'within_batch_chip',
    'within_chip_lane',
    'percent.mito'
)
covariate.columns <- c('cell', 'guide_count', 'prep_batch', 'percent.mito')
covariates <- covariates[, covariate.columns]

# read in cell cycle scores and merge with covariates
s.scores <- read.csv(
    'data/experimental/interim/cell_cycle_s_scores.csv'
)
g2m.scores <- read.csv(
    'data/experimental/interim/cell_cycle_g2m_scores.csv'
)
colnames(s.scores) <- c('cell', 's.score')
colnames(g2m.scores) <- c('cell', 'g2m.score')

covariates <- merge(covariates, s.scores)
covariates <- merge(covariates, g2m.scores)

# read in counts matrix
expr.matrix <- h5read(
    h5.name,
    'expr/expr_matrix'
)
genes <- h5read(
    h5.name,
    'expr/gene_names'
)
rownames(expr.matrix) <- genes
barcodes <- h5read(
    h5.name,
    'expr/cell_barcodes'
)
colnames(expr.matrix) <- barcodes

# compute scaling factors and add to covariates
scaling.factors <- data.frame(colSums(expr.matrix) / 1e6)
scaling.factors$cell <- rownames(scaling.factors)
rownames(scaling.factors) <- NULL
colnames(scaling.factors) <- c('scaling.factor', 'cell')
covariates <- merge(
    scaling.factors,
    covariates,
    by = 'cell',
    sort = FALSE
)

# write to h5 structure
h5write(
    covariates,
    h5.name,
    'expr/cell_covariates'
)

# read in guide matrix
guide.matrix <- h5read(
    'data/experimental/interim/guide_matrix.h5',
    'df/block0_values'
)
guide.matrix <- as.matrix(guide.matrix)
guide.matrix <- t(guide.matrix)

# define data dimensionality and chunk size
h5createDataset(
    h5.name,
    'grna/guide_matrix',
    c(nrow(guide.matrix), ncol(guide.matrix)),
    storage.mode = 'integer',
    chunk = c(2000, ncol(guide.matrix))
)

# create h5
h5write(
    guide.matrix,
    h5.name,
    'grna/guide_matrix'
)

# read in guide matrix guide names
guide.names <- h5read(
    'data/experimental/interim/guide_matrix.h5',
    'df/axis1'
)

# write guide names to h5 structure
h5write(
    guide.names,
    h5.name,
    'grna/guide_names'
)

# read in guide matrix cell barcodes
barcodes <- h5read(
    'data/experimental/interim/guide_matrix.h5',
    'df/axis0'
)

h5write(
    barcodes,
    h5.name,
    'grna/cell_barcodes'
)
