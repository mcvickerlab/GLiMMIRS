# This script runs the GLiMMIRS model on the 330 enhancer-enhancer pairs
# which were computed from the 664 previously published enhancer-gene pairs
# from Gasperini et al.
#
# Author: Karthik Guruvayurappan

library(rhdf5)
library(MASS)

# define h5 file name as a variable
h5.name <- 'data/experimental/processed/gasperini_data.h5'

# read in enhancer-enhancer pairs
enhancer.pairs <- h5read(
    h5.name,
    'enhancer_enhancer_330'
)



# read in covariates
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

# add pseudocount to count data
pseudocount <- 0.01
counts.matrix <- counts.matrix + pseudocount

# compute scaling factors based on count matrix
print('computing scaling factors!')
scaling.factors <- colSums(counts.matrix) / 1e6

# read in enhancer-enhancer pairs
enhancer.enhancer.pairs <- read.csv('/iblm/netapp/data1/external/Gasperini2019/processed/enhancer_pairs_suppl_table_2.csv')

intercept.coeff.list <- rep(NA, nrow(enhancer.enhancer.pairs))
intercept.pvalue.list <- rep(NA, nrow(enhancer.enhancer.pairs))
enhancer.1.list <- rep(NA, nrow(enhancer.enhancer.pairs))
enhancer.2.list <- rep(NA, nrow(enhancer.enhancer.pairs))
gene.list <- rep(NA, nrow(enhancer.enhancer.pairs))
enhancer.1.coeff.list <- rep(NA, nrow(enhancer.enhancer.pairs))
enhancer.1.pvalue.list <- rep(NA, nrow(enhancer.enhancer.pairs))
enhancer.2.coeff.list <- rep(NA, nrow(enhancer.enhancer.pairs))
enhancer.2.pvalue.list <- rep(NA, nrow(enhancer.enhancer.pairs))
interaction.coeff.list <- rep(NA, nrow(enhancer.enhancer.pairs))
interaction.pvalue.list <- rep(NA, nrow(enhancer.enhancer.pairs))
prep.batch.coeff.list <- rep(NA, nrow(enhancer.enhancer.pairs))
prep.batch.pvalue.list <- rep(NA, nrow(enhancer.enhancer.pairs))
guide.count.coeff.list <- rep(NA, nrow(enhancer.enhancer.pairs))
guide.count.pvalue.list <- rep(NA, nrow(enhancer.enhancer.pairs))
percent.mito.coeff.list <- rep(NA, nrow(enhancer.enhancer.pairs))
percent.mito.pvalue.list <- rep(NA, nrow(enhancer.enhancer.pairs))
s.score.coeff.list <- rep(NA, nrow(enhancer.enhancer.pairs))
s.score.pvalue.list <- rep(NA, nrow(enhancer.enhancer.pairs))
g2m.score.coeff.list <- rep(NA, nrow(enhancer.enhancer.pairs))
g2m.score.pvalue.list <- rep(NA, nrow(enhancer.enhancer.pairs))

for (i in 1:nrow(enhancer.enhancer.pairs)) {

    # get name of enhancers and gene
    enhancer.1 <- enhancer.enhancer.pairs[i, 'enhancer_1']
    enhancer.2 <- enhancer.enhancer.pairs[i, 'enhancer_2']
    gene <- enhancer.enhancer.pairs[i, 'gene']

    print(paste0('running model for ', enhancer.1, ' enhancer 1 and ', enhancer.2, ' enhancer 2 and ', gene, ' gene!'))

    enhancer.1.spacers <- enhancer.to.spacer.table[enhancer.to.spacer.table$target.site == enhancer.1, ]$spacer.sequence
    enhancer.2.spacers <- enhancer.to.spacer.table[enhancer.to.spacer.table$target.site == enhancer.2, ]$spacer.sequence

    # get guide effiencies corresponding to spacers
    enhancer.1.spacers.efficiencies <- guide.efficiencies.table[guide.efficiencies.table$spacer %in% enhancer.1.spacers, c('spacer', 'Cutting.Efficiency')]
    enhancer.2.spacers.efficiencies <- guide.efficiencies.table[guide.efficiencies.table$spacer %in% enhancer.2.spacers, c('spacer', 'Cutting.Efficiency')]

    enhancer.1.spacers.efficiencies[is.na(enhancer.1.spacers.efficiencies)] <- 0
    enhancer.2.spacers.efficiencies[is.na(enhancer.2.spacers.efficiencies)] <- 0

    enhancer.1.indicator.probs <- rep(1, nrow(cell.guide.matrix))

    for (j in 1:nrow(enhancer.1.spacers.efficiencies)) {
        guide.spacer <- enhancer.1.spacers.efficiencies$spacer[j]
        guide.efficiency <- enhancer.1.spacers.efficiencies$Cutting.Efficiency[j]
        guide.indicator.vector <- cell.guide.matrix[, guide.spacer]
        guide.probs <- 1 - (guide.indicator.vector * guide.efficiency)
        enhancer.1.indicator.probs <- enhancer.1.indicator.probs * guide.probs
    }

    enhancer.1.indicator.probs <- 1 - enhancer.1.indicator.probs
    enhancer.1.indicator.vector <- enhancer.1.indicator.probs

    enhancer.2.indicator.probs <- rep(1, nrow(cell.guide.matrix))

    for (j in 1:nrow(enhancer.2.spacers.efficiencies)) {
        guide.spacer <- enhancer.2.spacers.efficiencies$spacer[j]
        guide.efficiency <- enhancer.2.spacers.efficiencies$Cutting.Efficiency[j]
        guide.indicator.vector <- cell.guide.matrix[, guide.spacer]
        guide.probs <- 1 - (guide.indicator.vector * guide.efficiency)
        enhancer.2.indicator.probs <- enhancer.2.indicator.probs * guide.probs
    }

    enhancer.2.indicator.probs <- 1 - enhancer.2.indicator.probs
    enhancer.2.indicator.vector <- enhancer.2.indicator.probs

    # get gene counts for gene
    gene.counts <- counts.matrix[gene, ]

    # create dataframe for modeling
    model.df <- cbind(covariates, enhancer.1.indicator.vector, enhancer.2.indicator.vector, gene.counts)

    # fit negative binomial GLM model
    model <- glm.nb(
        formula = gene.counts ~ enhancer.1.indicator.vector * enhancer.2.indicator.vector + prep_batch + guide_count + percent.mito + s.score + g2m.score + offset(log(scaling.factors)),
        data = model.df
    )

    enhancer.1.list[i] <- enhancer.1
    enhancer.2.list[i] <- enhancer.2
    gene.list[i] <- gene
    if ('(Intercept)' %in% rownames(summary(model)$coefficients)){
        intercept.coeff.list[i] <- summary(model)$coefficients['(Intercept)', 'Estimate']
        intercept.pvalue.list[i] <- summary(model)$coefficients['(Intercept)', 'Pr(>|z|)']

    }
    else {
        intercept.coeff.list[i] <- NA
        intercept.pvalue.list[i] <- NA
    }
    if ('enhancer.1.indicator.vector' %in% rownames(summary(model)$coefficients)){
        enhancer.1.coeff.list[i] <- summary(model)$coefficients['enhancer.1.indicator.vector', 'Estimate']
        enhancer.1.pvalue.list[i] <- summary(model)$coefficients['enhancer.1.indicator.vector', 'Pr(>|z|)']

    }
    else {
        enhancer.1.coeff.list[i] <- NA
        enhancer.1.pvalue.list[i] <- NA
    }

    if ('enhancer.2.indicator.vector' %in% rownames(summary(model)$coefficients)){
        enhancer.2.coeff.list[i] <- summary(model)$coefficients['enhancer.2.indicator.vector', 'Estimate']
        enhancer.2.pvalue.list[i] <- summary(model)$coefficients['enhancer.2.indicator.vector', 'Pr(>|z|)']

    }
    else {
        enhancer.2.coeff.list[i] <- NA
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

    if ('prep_batch' %in% rownames(summary(model)$coefficients)){
        prep.batch.coeff.list[i] <- summary(model)$coefficients['prep_batch', 'Estimate']
        prep.batch.pvalue.list[i] <- summary(model)$coefficients['prep_batch', 'Pr(>|z|)']

    }
    else {
        prep.batch.coeff.list[i] <- NA
        prep.batch.pvalue.list[i] <- NA
    }

    if ('guide_count' %in% rownames(summary(model)$coefficients)){
        guide.count.coeff.list[i] <- summary(model)$coefficients['guide_count', 'Estimate']
        guide.count.pvalue.list[i] <- summary(model)$coefficients['guide_count', 'Pr(>|z|)']

    }
    else {
        guide.count.coeff.list[i] <- NA
        guide.count.pvalue.list[i] <- NA
    }

    if ('percent.mito' %in% rownames(summary(model)$coefficients)){
        percent.mito.coeff.list[i] <- summary(model)$coefficients['percent.mito', 'Estimate']
        percent.mito.pvalue.list[i] <- summary(model)$coefficients['percent.mito', 'Pr(>|z|)']

    }
    else {
        percent.mito.coeff.list[i] <- NA
        percent.mito.pvalue.list[i] <- NA
    }

    if ('s.score' %in% rownames(summary(model)$coefficients)){
        s.score.coeff.list[i] <- summary(model)$coefficients['s.score', 'Estimate']
        s.score.pvalue.list[i] <- summary(model)$coefficients['s.score', 'Pr(>|z|)']

    }
    else {
        s.score.coeff.list[i] <- NA
        s.score.pvalue.list[i] <- NA
    }

    if ('g2m.score' %in% rownames(summary(model)$coefficients)){
        g2m.score.coeff.list[i] <- summary(model)$coefficients['g2m.score', 'Estimate']
        g2m.score.pvalue.list[i] <- summary(model)$coefficients['g2m.score', 'Pr(>|z|)']

    }
    else {
        g2m.score.coeff.list[i] <- NA
        g2m.score.pvalue.list[i] <- NA
    }
}

# write to output file
print('writing p-values to output file!')
pvalue.table <- cbind(
    intercept.coeff.list, intercept.pvalue.list,
    enhancer.1.list, enhancer.2.list, gene.list,
    enhancer.1.coeff.list, enhancer.1.pvalue.list,
    enhancer.2.coeff.list, 
    enhancer.2.pvalue.list,
    interaction.coeff.list, interaction.pvalue.list, 
    prep.batch.coeff.list, prep.batch.pvalue.list,
    guide.count.coeff.list, guide.count.pvalue.list,
    percent.mito.coeff.list, percent.mito.pvalue.list,
    s.score.coeff.list, s.score.pvalue.list,
    g2m.score.coeff.list, g2m.score.pvalue.list)

write.csv(
    pvalue.table,
    '/iblm/netapp/data1/external/Gasperini2019/processed/23_10_31_enhancer_enhancer_pairs_suppl_table_2_pseudocount_model_enhancer_effects.csv',
    row.names = FALSE
)
