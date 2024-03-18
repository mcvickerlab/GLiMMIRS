# This script splits the at-scale enhancer pairs into 32 batches so that
# GLiMMIRS can be parallelized.
#
# Author: Karthik Guruvayurappan

library(rhdf5)

# read in at-scale enhancer pairs
h5.name <- 'data/experimental/processed/gasperini_data.h5'

enhancer.pairs <- read.csv(
    'data/experimental/processed/perturbation_counts.csv'
)

# read in gene names (to filter pairs)
genes <- h5read(
    h5.name,
    'expr/gene_names'
)
enhancer.pairs <- enhancer.pairs[enhancer.pairs$gene %in% genes, ]

# filter to >= 10 cell threshold
enhancer.pairs <- enhancer.pairs[enhancer.pairs$double.count.list >= 10, ]

enhancer.pairs <- enhancer.pairs[, c('enh.1.list', 'enh.2.list', 'gene.list')]
colnames(enhancer.pairs) <- c('enhancer_1', 'enhancer_2', 'gene')

# dimensions should be 82,314 x 6 for the enhancer pairs
rownames(enhancer.pairs) <- NULL
batches <- as.integer(rownames(enhancer.pairs)) %% 32
enhancer.pairs <- split(enhancer.pairs, batches)

# write out batch enhancer pair files
for (i in 1:32) {
    write.csv(
        enhancer.pairs[[i]],
        paste0('data/experimental/interim/enhancer_pairs_', i, '.csv'),
        row.names = FALSE,
        quote = FALSE
    )
}
