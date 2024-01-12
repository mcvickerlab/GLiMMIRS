# This script runs our "baseline" enhancer-gene model on the 664 enhancer-gene
# pairs previously published by Gasperini et al. Our models differ by including
# guide efficiency information and cell cycle scores.
#
# Author: Karthik Guruvayurappan

library(MASS)
library(rhdf5)

# define h5 file name as a variable
h5.name <- 'data/experimental/processed/gasperini_data.h5'

# read in enhancer-gene pairs
enhancer.gene <- h5read(
    h5.name,
    'enhancer_gene'
)

# read in enhancer-guide table
enhancer.guide <- h5read(
    h5.name,
    'enhancer_guide'
)

# read in guide efficiency information
guide.info <- h5read(
    h5.name,
    'grna/guide_info'
)

# read in guide matrix
guide.matrix <- h5read(
    h5.name,
    'grna/guide_matrix'
)
guide.names <- h5read(
    h5.name,
    'grna/guide_names'
)
rownames(guide.matrix) <- guide.names
barcodes <- h5read(
    h5.name,
    'grna/cell_barcodes'
)
colnames(guide.matrix) <- barcodes

# read in cell covariates
covariates <- h5read(
    h5.name,
    'expr/cell_covariates'
)

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

# add pseudocount to count data
pseudocount <- 0.01
expr.matrix <- expr.matrix + pseudocount

# compute scaling factors based on count matrix
# and add to covariates
# the cell barcodes for the covariates and expression matrix will not match
# otherwise!
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








enhancer.list <- rep(NA, nrow(enhancer.gene.pairs))
gene.list <- rep(NA, nrow(enhancer.gene.pairs))
pvalue.list <- rep(NA, nrow(enhancer.gene.pairs))
enhancer.effect.list <- rep(NA, nrow(enhancer.gene.pairs))

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
    gene.counts <- counts.matrix[gene, ]

    # create dataframe for modeling
    model.df <- cbind(covariates, indicator.vector, gene.counts)

    # fit negative binomial GLM mode
    model <- glm.nb(
        formula = gene.counts ~ indicator.vector + prep_batch + guide_count + percent.mito + s.score + g2m.score + offset(log(scaling.factors)),
        data = model.df
    )

    enhancer.list[i] <- enhancer
    gene.list[i] <- gene

    if ('indicator.vector' %in% rownames(summary(model)$coefficients)){
        enhancer.effect.list[i] <- summary(model)$coefficients['indicator.vector', 'Estimate']
        pvalue.list[i] <- summary(model)$coefficients['indicator.vector', 'Pr(>|z|)']

    }
    else {
        enhancer.effect.list[i] <- NA
        pvalue.list[i] <- NA
    }
}

# write to output file
print('writing p-values to output file!')
pvalue.table <- cbind(enhancer.list, gene.list, enhancer.effect.list, pvalue.list)
write.csv(
    pvalue.table,
    '/iblm/netapp/data1/external/Gasperini2019/processed/23_11_19_enhancer_gene_pairs_suppl_table_2_baseline_pseudocount_model.csv',
    row.names = FALSE
)
