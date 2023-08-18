# This program computes the number of cells which contains guides for both enhancers in the
# at-scale enhancer pairs determined from the Gasperini et al. dataset (roughly 1M)
#
# Author: Karthik Guruvayurappan

library(rhdf5)

# read in table mapping enhancers to spacers and reformat enhancer names
print('reading in enhancer-to-spacer table!')
enhancer.to.spacer.table <- read.table(
    '/iblm/netapp/data1/external/Gasperini2019/suppl/GSE120861_grna_groups.at_scale.txt',
    sep = '\t'
)
colnames(enhancer.to.spacer.table) <- c('target.site', 'spacer.sequence')
enhancer.to.spacer.table$target.site <- sapply(enhancer.to.spacer.table$target.site, FUN = function(x) {
    if (startsWith(x, 'chr')) {
        return (strsplit(x, '_')[[1]][1])
    }
    else {
        return (x)
    }
})

# read in cell-guide matrix
print('reading in cell-guide matrix!')
cell.guide.matrix <- h5read('/iblm/netapp/data1/external/Gasperini2019/processed/gasperini_data.h5', 'cell.guide.matrix')
guide.spacers <- h5read('/iblm/netapp/data1/external/Gasperini2019/processed/gasperini_data.h5', 'guide.spacers')
colnames(cell.guide.matrix) <- guide.spacers

enhancer.enhancer.pairs <- read.csv('/iblm/netapp/data1/external/Gasperini2019/processed/at_scale_enhancer_enhancer_pairs.csv')
enhancer.enhancer.pairs$enhancer_1 <- sapply(enhancer.enhancer.pairs$enhancer_1, FUN = function(x) {
    if (startsWith(x, 'chr')) {
        return (strsplit(x, '_')[[1]][1])
    }
    else {
        return (x)
    }
})
enhancer.enhancer.pairs$enhancer_2 <- sapply(enhancer.enhancer.pairs$enhancer_2, FUN = function(x) {
    if (startsWith(x, 'chr')) {
        return (strsplit(x, '_')[[1]][1])
    }
    else {
        return (x)
    }
})

# drop duplicates from enhancer pairs
enhancer.enhancer.pairs <- enhancer.enhancer.pairs[!duplicated(enhancer.enhancer.pairs), ]

# filter for enhancer-enhancer pairs with gene expression info
gene.names <- h5read('/iblm/netapp/data1/external/Gasperini2019/processed/gasperini_data.h5', 'gene.names')
enhancer.enhancer.pairs <- enhancer.enhancer.pairs[enhancer.enhancer.pairs$gene %in% gene.names, ]

enhancer.1.list <- rep(NA, nrow(enhancer.enhancer.pairs))
enhancer.2.list <- rep(NA, nrow(enhancer.enhancer.pairs))
enhancer.1.count.list <- rep(NA, nrow(enhancer.enhancer.pairs))
enhancer.2.count.list <- rep(NA, nrow(enhancer.enhancer.pairs))
count.list <- rep(NA, nrow(enhancer.enhancer.pairs))

for (i in 1:nrow(enhancer.enhancer.pairs)) {

    # get name of enhancers and gene
    enhancer.1 <- enhancer.enhancer.pairs[i, 'enhancer_1']
    enhancer.2 <- enhancer.enhancer.pairs[i, 'enhancer_2']

    enhancer.1.spacers <- enhancer.to.spacer.table[enhancer.to.spacer.table$target.site == enhancer.1, ]$spacer.sequence
    enhancer.2.spacers <- enhancer.to.spacer.table[enhancer.to.spacer.table$target.site == enhancer.2, ]$spacer.sequence

    enhancer.1.spacers <- enhancer.1.spacers[enhancer.1.spacers %in% colnames(cell.guide.matrix)]
    enhancer.2.spacers <- enhancer.2.spacers[enhancer.2.spacers %in% colnames(cell.guide.matrix)]

    enhancer.1.list[i] <- enhancer.1
    enhancer.2.list[i] <- enhancer.2

    enhancer.1.indicator.vector <- rep(0, nrow(cell.guide.matrix))

    for (j in 1:length(enhancer.1.spacers)) {
        guide.spacer <- enhancer.1.spacers[j]
        guide.indicator.vector <- cell.guide.matrix[, guide.spacer]
        enhancer.1.indicator.vector <- enhancer.1.indicator.vector + guide.indicator.vector
    }

    enhancer.1.indicator.vector[enhancer.1.indicator.vector > 1] <- 1

    enhancer.2.indicator.vector <- rep(0, nrow(cell.guide.matrix))

    for (j in 1:length(enhancer.2.spacers)) {
        guide.spacer <- enhancer.2.spacers[j]
        guide.indicator.vector <- cell.guide.matrix[, guide.spacer]
        enhancer.2.indicator.vector <- enhancer.2.indicator.vector + guide.indicator.vector
    }

    enhancer.2.indicator.vector[enhancer.2.indicator.vector > 1] <- 1

    sum.vector <- enhancer.1.indicator.vector + enhancer.2.indicator.vector

    both.target.count <- sum(sum.vector == 2)
    count.list[i] <- both.target.count

    enhancer.1.count.list[i] <- sum(enhancer.1.indicator.vector) - both.target.count
    enhancer.2.count.list[i] <- sum(enhancer.2.indicator.vector) - both.target.count
}

# write to output file
print('writing counts to output file!')
count.table <- cbind(enhancer.1.list, enhancer.2.list, enhancer.1.count.list, enhancer.2.count.list, count.list)
write.csv(
    count.table,
    '/iblm/netapp/data1/external/Gasperini2019/processed/23_08_03_at_scale_enhancer_enhancer_pairs_cells_counts.csv',
    row.names = FALSE
)
