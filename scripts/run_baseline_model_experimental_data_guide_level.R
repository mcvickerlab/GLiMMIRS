# This program runs the baseline GLM model on all 664 enhancer-gene pairs previously published by
# Gasperini et al. in 2019. However, this program separates the indicator vectors by guide, with
# the intention of demonstrating why guide efficiency modeling is necessary. This program was 
# written by Karthik Guruvayurappan.

library(MASS)
library(rhdf5)

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

# read in counts matrix
print('reading in counts matrix!')
counts.matrix <- h5read('/iblm/netapp/data1/external/Gasperini2019/processed/gasperini_data.h5', 'gene.counts')
gene.names <- h5read('/iblm/netapp/data1/external/Gasperini2019/processed/gasperini_data.h5', 'gene.names')
rownames(counts.matrix) <- gene.names

# compute scaling factors based on count matrix
print('computing scaling factors!')
scaling.factors <- colSums(counts.matrix) / 1e6

# read in previously published enhancer gene pairs
print('reading in enhancer gene pairs!')
enhancer.gene.pairs <- read.csv(
    '/iblm/netapp/data1/external/Gasperini2019/gasperini_enhancer_gene_pairs_suppl_table_2.csv'
)

enhancer.list <- rep(NA, nrow(enhancer.gene.pairs))
gene.list <- rep(NA, nrow(enhancer.gene.pairs))
efficiency.range.list <- rep(NA, nrow(enhancer.gene.pairs))
effect.range.list <- rep(NA, nrow(enhancer.gene.pairs))

for (i in 1:nrow(enhancer.gene.pairs)) {

    # get name of enhancer and gene
    enhancer <- enhancer.gene.pairs[i, 'Target_Site']
    gene <- enhancer.gene.pairs[i, 'ENSG']

    print(paste0('running model for ', enhancer, ' enhancer and ', gene, ' gene!'))

    # get spacer sequences corresponding to enhancer
    enhancer.spacers <- enhancer.to.spacer.table[enhancer.to.spacer.table$target.site == enhancer, ]$spacer.sequence
    
    # get guide effiencies corresponding to spacers
    enhancer.spacers.efficiencies <- guide.efficiencies.table[guide.efficiencies.table$spacer %in% enhancer.spacers, c('spacer', 'Cutting.Efficiency')]
    enhancer.spacers.efficiencies[is.na(enhancer.spacers.efficiencies)] <- 0

    guide.indicator.vectors <- cell.guide.matrix[, enhancer.spacers.efficiencies$spacer]

    # get gene counts for gene
    gene.counts <- counts.matrix[gene, ]

    # create dataframe for modeling
    model.df <- cbind(covariates, guide.indicator.vectors, gene.counts)

    # fit negative binomial GLM model
    guide.formula = 'gene.counts ~ '
    for (j in 1:nrow(enhancer.spacers.efficiencies)) {
        guide.formula <- paste0(guide.formula, enhancer.spacers.efficiencies$spacer[j], ' + ')
    }

    guide.formula <- paste0(guide.formula, 'prep_batch + guide_count + percent.mito + s.score + g2m.score + offset(scaling.factors)')

    model <- glm.nb(
        formula = as.formula(guide.formula),
        data = model.df
    )
    
    enhancer.list[i] <- enhancer
    gene.list[i] <- gene

    efficiency.range <- max(enhancer.spacers.efficiencies$Cutting.Efficiency) - min(enhancer.spacers.efficiencies$Cutting.Efficiency)
    efficiency.range.list[i] <- efficiency.range

    coefficient.stats <- summary(model)$coefficients[enhancer.spacers.efficiencies$spacer, ]
    effect.range <- max(coefficient.stats[, 'Estimate']) - min(coefficient.stats[, 'Estimate'])
    effect.range.list[i] <- effect.range
}

# write to output file
print('writing p-values to output file!')
effect.table <- cbind(enhancer.list, gene.list, efficiency.range.list, effect.range.list)
write.csv(
    effect.table,
    '/iblm/netapp/data1/external/Gasperini2019/processed/enhancer_gene_pairs_suppl_table_2_baseline_model_guide_level.csv',
    row.names = FALSE
)
