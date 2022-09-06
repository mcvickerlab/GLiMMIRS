# This program computes the number of cells which would have both enhancers targeted for the 330
# enhancer-enhancer pairs derived from the 664 enhancer-gene pairs published by Gasperini et al.
# in 2019.
#
# Author: Karthik Guruvayurappan

library(rhdf5)
library(MASS)

# read in and sort covariates
print('reading in covariates!')
covariates <- h5read(
    file = '/iblm/netapp/data1/external/Gasperini2019/processed/gasperini_data.h5',
    name = 'covariates'
)
cell.barcodes <- h5read(
    file = '/iblm/netapp/data1/external/Gasperini2019/processed/gasperini_data.h5',
    name = 'cell.barcodes'
)
covariates <- merge(
    data.frame(cell.barcodes),
    covariates,
    by.x = 'cell.barcodes',
    by.y = 'cell',
    sort = FALSE
)

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

# read in guide efficiency information
print('reading in guide efficiencies!')
guide.efficiencies.table <- h5read(
    '/iblm/netapp/data1/external/Gasperini2019/processed/gasperini_data.h5',
    'guidescan.output'
)
guide.efficiencies.table$spacer <- substring(
    guide.efficiencies.table$gRNA,
    1,
    nchar(guide.efficiencies.table$gRNA) - 3
)

# read in cell-guide matrix
print('reading in cell-guide matrix!')
cell.guide.matrix <- h5read('/iblm/netapp/data1/external/Gasperini2019/processed/gasperini_data.h5', 'cell.guide.matrix')
guide.spacers <- h5read('/iblm/netapp/data1/external/Gasperini2019/processed/gasperini_data.h5', 'guide.spacers')
colnames(cell.guide.matrix) <- guide.spacers

# read in enhancer-enhancer pairs
enhancer.enhancer.pairs <- read.csv('/iblm/netapp/data1/external/Gasperini2019/processed/enhancer_pairs_suppl_table_2.csv')

enhancer.1.list <- rep(NA, nrow(enhancer.enhancer.pairs))
enhancer.2.list <- rep(NA, nrow(enhancer.enhancer.pairs))
count.list <- rep(NA, nrow(enhancer.enhancer.pairs))
count.guide.efficiency.list <- rep(NA, nrow(enhancer.enhancer.pairs))

for (i in 1:nrow(enhancer.enhancer.pairs)) {

    # get name of enhancers and gene
    enhancer.1 <- enhancer.enhancer.pairs[i, 'enhancer_1']
    enhancer.2 <- enhancer.enhancer.pairs[i, 'enhancer_2']

    enhancer.1.spacers <- enhancer.to.spacer.table[enhancer.to.spacer.table$target.site == enhancer.1, ]$spacer.sequence
    enhancer.2.spacers <- enhancer.to.spacer.table[enhancer.to.spacer.table$target.site == enhancer.2, ]$spacer.sequence

    # get guide effiencies corresponding to spacers
    enhancer.1.spacers.efficiencies <- guide.efficiencies.table[guide.efficiencies.table$spacer %in% enhancer.1.spacers, c('spacer', 'Cutting.Efficiency')]
    enhancer.2.spacers.efficiencies <- guide.efficiencies.table[guide.efficiencies.table$spacer %in% enhancer.2.spacers, c('spacer', 'Cutting.Efficiency')]

    enhancer.1.spacers.efficiencies[is.na(enhancer.1.spacers.efficiencies)] <- 0
    enhancer.2.spacers.efficiencies[is.na(enhancer.2.spacers.efficiencies)] <- 0

    enhancer.1.indicator.vector <- rep(0, nrow(cell.guide.matrix))

    for (j in 1:nrow(enhancer.1.spacers.efficiencies)) {
        guide.spacer <- enhancer.1.spacers.efficiencies$spacer[j]
        guide.indicator.vector <- cell.guide.matrix[, guide.spacer]

        enhancer.1.indicator.vector <- enhancer.1.indicator.vector + guide.indicator.vector
    }

    enhancer.1.indicator.vector[enhancer.1.indicator.vector > 1] <- 1 

    enhancer.2.indicator.vector <- rep(0, nrow(cell.guide.matrix))

    for (j in 1:nrow(enhancer.2.spacers.efficiencies)) {
        guide.spacer <- enhancer.2.spacers.efficiencies$spacer[j]
        guide.indicator.vector <- cell.guide.matrix[, guide.spacer]
        
        enhancer.2.indicator.vector <- enhancer.2.indicator.vector + guide.indicator.vector
    }

    enhancer.2.indicator.vector[enhancer.2.indicator.vector > 1] <- 1 
    
    both.target.count <- sum((enhancer.1.indicator.vector + enhancer.2.indicator.vector) == 2)

    enhancer.1.list[i] <- enhancer.1
    enhancer.2.list[i] <- enhancer.2
    count.list[i] <- both.target.count

    enhancer.1.indicator.probs <- rep(1, nrow(cell.guide.matrix))

    for (j in 1:nrow(enhancer.1.spacers.efficiencies)) {
        guide.spacer <- enhancer.1.spacers.efficiencies$spacer[j]
        guide.efficiency <- enhancer.1.spacers.efficiencies$Cutting.Efficiency[j]
        guide.indicator.vector <- cell.guide.matrix[, guide.spacer]
        guide.probs <- 1 - (guide.indicator.vector * guide.efficiency)
        enhancer.1.indicator.probs <- enhancer.1.indicator.probs * guide.probs
    }

    enhancer.1.indicator.probs <- 1 - enhancer.1.indicator.probs
    enhancer.1.indicator.vector <- rbinom(nrow(cell.guide.matrix), 1, enhancer.1.indicator.probs)

    enhancer.2.indicator.probs <- rep(1, nrow(cell.guide.matrix))

    for (j in 1:nrow(enhancer.2.spacers.efficiencies)) {
        guide.spacer <- enhancer.2.spacers.efficiencies$spacer[j]
        guide.efficiency <- enhancer.2.spacers.efficiencies$Cutting.Efficiency[j]
        guide.indicator.vector <- cell.guide.matrix[, guide.spacer]
        guide.probs <- 1 - (guide.indicator.vector * guide.efficiency)
        enhancer.2.indicator.probs <- enhancer.2.indicator.probs * guide.probs
    }

    enhancer.2.indicator.probs <- 1 - enhancer.2.indicator.probs
    enhancer.2.indicator.vector <- rbinom(nrow(cell.guide.matrix), 1, enhancer.2.indicator.probs)

    both.target.count <- sum((enhancer.1.indicator.vector + enhancer.2.indicator.vector) == 2)
    count.guide.efficiency.list[i] <- both.target.count


}

# write to output file
print('writing counts to output file!')
pvalue.table <- cbind(enhancer.1.list, enhancer.2.list, count.list, count.guide.efficiency.list)
write.csv(
    pvalue.table,
    '/iblm/netapp/data1/external/Gasperini2019/processed/enhancer_enhancer_pairs_suppl_table_2_count_both_enhancers_target.csv',
    row.names = FALSE
)
