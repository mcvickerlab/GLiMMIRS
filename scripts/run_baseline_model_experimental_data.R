# This program runs the baseline GLM model on all 664 enhancer-gene pairs previously published by
# Gasperini et al. in 2019. This program was written by Karthik Guruvayurappan.

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
pvalue.list <- rep(NA, nrow(enhancer.gene.pairs))

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

    indicator.vector.probs <- rep(1, nrow(cell.guide.matrix))

    for (j in 1:nrow(enhancer.spacers.efficiencies)) {
        
        guide.spacer <- enhancer.spacers.efficiencies$spacer[j]
        guide.efficiency <- enhancer.spacers.efficiencies$Cutting.Efficiency[j]

        guide.indicator.vector <- cell.guide.matrix[, guide.spacer]

        guide.probs <- 1 - (guide.indicator.vector * guide.efficiency)

        indicator.vector.probs <- indicator.vector.probs * guide.probs
    }

    indicator.vector.probs <- 1 - indicator.vector.probs
    indicator.vector <- indicator.vector.probs

    # get gene counts for gene
    pseudocount <- 0.01
    gene.counts <- counts.matrix[gene, ] + pseudocount

    # create dataframe for modeling
    model.df <- cbind(covariates, indicator.vector, gene.counts)

    # fit negative binomial GLM model
    model <- glm.nb(
        formula = gene.counts ~ indicator.vector + prep_batch + guide_count + percent.mito + s.score + g2m.score + offset(scaling.factors),
        data = model.df
    )

    enhancer.list[i] <- enhancer
    gene.list[i] <- gene

    if ('indicator.vector' %in% rownames(summary(model)$coefficients)){
        pvalue.list[i] <- summary(model)$coefficients['indicator.vector', 'Pr(>|z|)']

    }
    else {
        pvalue.list[i] <- NA
    }
}

# write to output file
print('writing p-values to output file!')
pvalue.table <- cbind(enhancer.list, gene.list, pvalue.list)
write.csv(
    pvalue.table,
    '/iblm/netapp/data1/external/Gasperini2019/processed/enhancer_gene_pairs_suppl_table_2_baseline_pseudocount_model.csv',
    row.names = FALSE
)
