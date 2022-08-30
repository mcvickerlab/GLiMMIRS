# This program runs a generalized linear model on enhancer pairs determined by analyzing
# the 664 enhancer-gene pairs published in the Gasperini et al. 2019 paper, and looking at
# enhancers that target the same gene.
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

# read in counts matrix
print('reading in counts matrix!')
counts.matrix <- h5read('/iblm/netapp/data1/external/Gasperini2019/processed/gasperini_data.h5', 'gene.counts')
gene.names <- h5read('/iblm/netapp/data1/external/Gasperini2019/processed/gasperini_data.h5', 'gene.names')
rownames(counts.matrix) <- gene.names

# compute scaling factors based on count matrix
print('computing scaling factors!')
scaling.factors <- colSums(counts.matrix) / 1e6

# read in enhancer-enhancer pairs
enhancer.enhancer.pairs <- read.csv('/iblm/netapp/data1/external/Gasperini2019/processed/enhancer_pairs_suppl_table_2.csv')

enhancer.1.list <- rep(NA, nrow(enhancer.enhancer.pairs))
enhancer.2.list <- rep(NA, nrow(enhancer.enhancer.pairs))
gene.list <- rep(NA, nrow(enhancer.enhancer.pairs))
enhancer.1.pvalue.list <- rep(NA, nrow(enhancer.enhancer.pairs))
enhancer.2.pvalue.list <- rep(NA, nrow(enhancer.enhancer.pairs))
interaction.coeff.list <- rep(NA, nrow(enhancer.enhancer.pairs))
interaction.pvalue.list <- rep(NA, nrow(enhancer.enhancer.pairs))

for (i in 1:nrow(enhancer.enhancer.pairs)) {

    # get name of enhancers and gene
    enhancer.1 <- enhancer.enhancer.pairs[i, 'enhancer_1']
    enhancer.2 <- enhancer.enhancer.pairs[i, 'enhancer_2']
    gene <- enhancer.enhancer.pairs[i, 'gene']

    print(paste0('running model for ', enhancer.1, ' enhancer 1 and ', enhancer.2, ' enhancer 2 and ', gene, ' gene!'))

    enhancer.1.spacers <- target.site.spacers.table[target.site.spacers.table$target.site == enhancer.1, ]$spacer.sequence
    enhancer.2.spacers <- target.site.spacers.table[target.site.spacers.table$target.site == enhancer.2, ]$spacer.sequence

    # get guide effiencies corresponding to spacers
    enhancer.1.spacers.efficiencies <- all.guide.efficiencies[all.guide.efficiencies$spacer %in% enhancer.1.spacers, c('spacer', 'Cutting.Efficiency')]
    enhancer.2.spacers.efficiencies <- all.guide.efficiencies[all.guide.efficiencies$spacer %in% enhancer.2.spacers, c('spacer', 'Cutting.Efficiency')]

    enhancer.1.spacers.efficiencies[is.na(enhancer.1.spacers.efficiencies)] <- 0
    enhancer.2.spacers.efficiencies[is.na(enhancer.2.spacers.efficiencies)] <- 0

    enhancer.1.indicator.probs <- rep(1, nrow(cell.guide.matrix))

    for (i in 1:nrow(enhancer.1.spacers.efficiencies)) {
        guide.spacer <- enhancer.1.spacers.efficiencies$spacer[i]
        guide.efficiency <- enhancer.1.spacers.efficiencies$Cutting.Efficiency[i]
        guide.indicator.vector <- cell.guide.matrix[, guide.spacer]
        guide.probs <- 1 - (guide.indicator.vector * guide.efficiency)
        enhancer.1.indicator.probs <- enhancer.1.indicator.probs * guide.probs
    }

    enhancer.1.indicator.probs <- 1 - enhancer.1.indicator.probs
    enhancer.1.indicator.vector <- rbinom(nrow(cell.guide.matrix), 1, enhancer.1.indicator.probs)

    enhancer.2.indicator.probs <- rep(1, nrow(cell.guide.matrix))

    for (i in 1:nrow(enhancer.2.spacers.efficiencies)) {
        guide.spacer <- enhancer.2.spacers.efficiencies$spacer[i]
        guide.efficiency <- enhancer.2.spacers.efficiencies$Cutting.Efficiency[i]
        guide.indicator.vector <- cell.guide.matrix[, guide.spacer]
        guide.probs <- 1 - (guide.indicator.vector * guide.efficiency)
        enhancer.2.indicator.probs <- enhancer.2.indicator.probs * guide.probs
    }

    enhancer.2.indicator.probs <- 1 - enhancer.2.indicator.probs
    enhancer.2.indicator.vector <- rbinom(nrow(cell.guide.matrix), 1, enhancer.2.indicator.probs)

    # get gene counts for gene
    gene.counts <- counts.matrix[gene, ]

    # create dataframe for modeling
    model.df <- cbind(covariates, enhancer.1.indicator.vector, enhancer.2.indicator.vector, gene.counts)

    # fit negative binomial GLM model
    model <- glm.nb(
        formula = gene.counts ~ enhancer.1.indicator.vector * enhancer.2.indicator.vector + prep_batch + guide_count + percent.mito + s.score + g2m.score + offset(scaling.factors),
        data = model.df
    )

    enhancer.1.list[i] <- enhancer.1
    enhancer.2.list[i] <- enhancer.2
    gene.list[i] <- gene

    if ('enhancer.1.indicator.vector' %in% rownames(summary(model)$coefficients)){
        enhancer.1.pvalue.list[i] <- summary(model)$coefficients['enhancer.1.indicator.vector', 'Pr(>|z|)']

    }
    else {
        enhancer.1.pvalue.list[i] <- NA
    }

    if ('enhancer.2.indicator.vector' %in% rownames(summary(model)$coefficients)){
        enhancer.2.pvalue.list[i] <- summary(model)$coefficients['enhancer.2.indicator.vector', 'Pr(>|z|)']

    }
    else {
        enhancer.2.pvalue.list[i] <- NA
    }

    if ('enhancer.1.indicator.vector:enhancer.2.indicator.vector' %in% rownames(summary(model)$coefficients)){
        interaction.coeff.list[i] <- summary(model)$coefficients['enhancer.1.indicator.vector:enhancer.2.indicator.vector', 'Estimate']
        interaction.pvalue.list[i] <- summary(model)$coefficients['enhancer.1.indicator.vector:enhancer.2.indicator.vector', 'Pr(>|z|)']

    }
    else {
        interaction.coeff.list[i] <- NA
        interaction.pvalue.list[i] <- NA
    }
}

# write to output file
print('writing p-values to output file!')
pvalue.table <- cbind(enhancer.1.list, enhancer.2.list, gene.list, enhancer.1.pvalue.list, enhancer.2.pvalue.list, interaction.coeff.list, interaction.pvalue.list)
write.csv(
    pvalue.table,
    '/iblm/netapp/data1/external/Gasperini2019/processed/enhancer_enhancer_pairs_suppl_table_2_model.csv',
    row.names = FALSE
)
