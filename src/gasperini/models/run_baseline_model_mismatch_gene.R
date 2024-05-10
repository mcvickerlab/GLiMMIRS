# This script runs our "baseline" enhancer-gene model on the 664 enhancer-gene
# pairs previously published by Gasperini et al. Our models differ by including
# guide efficiency information and cell cycle scores.
#
# Author: Karthik Guruvayurappan

library(MASS)
library(rhdf5)

set.seed(1)

# define h5 file name as a variable
h5.name <- 'data/gasperini/processed/gasperini_data.h5'

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
covariates$scaling.factor <- as.vector(covariates$scaling.factor)

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

enhancer.list <- rep(NA, nrow(enhancer.gene))
gene.list <- rep(NA, nrow(enhancer.gene))
effect.list <- rep(NA, nrow(enhancer.gene))
pvalue.list <- rep(NA, nrow(enhancer.gene))

for (i in 1:nrow(enhancer.gene)) {

    # get name of enhancer and gene
    enhancer <- enhancer.gene[i, 'Target_Site']
    gene <- sample(enhancer.gene$ENSG, 1)

    # print statement (for progress)
    print(paste0('running ', enhancer, ' and ', gene, '!'))

    # get enhancer spacer sequences
    spacers <- enhancer.guide[enhancer.guide$target.site == enhancer, 'spacer']

    # get spacer guide efficiencies and set NAs to 0 (excludes)
    efficiencies <- guide.info[guide.info$spacer %in% spacers, ]
    efficiencies <- efficiencies[, c('spacer', 'Cutting.Efficiency')]
    efficiencies[is.na(efficiencies)] <- 0

    # compute guide perturbation vector by computing 1 - (no perturbation prob)
    no.perturbation <- rep(1, ncol(guide.matrix))
    for (j in 1:nrow(efficiencies)) {
        spacer <- efficiencies[j, 'spacer']
        efficiency <- efficiencies[j, 'Cutting.Efficiency']
        spacer.vector <- guide.matrix[spacer, ]
        spacer.perturbation <- spacer.vector * efficiency
        spacer.no.perturbation <- 1 - spacer.perturbation
        no.perturbation <- no.perturbation * spacer.no.perturbation
    }
    perturbation <- 1 - no.perturbation

    # get gene counts
    gene.counts <- expr.matrix[gene, ]

    # fit negative binomial GLM
    model.df <- cbind(perturbation, gene.counts, covariates)
    model.formula <- as.formula(paste0(
        'gene.counts ~ ',
        'perturbation + ',
        'prep_batch + ',
        'guide_count + ',
        'percent.mito + ',
        's.score + ',
        'g2m.score + ',
        'offset(log(scaling.factor))'
    ))
    model <- glm.nb(
        formula = model.formula,
        data = model.df
    )

    # store model outputs
    enhancer.list[i] <- enhancer
    gene.list[i] <- gene
    if ('perturbation' %in% rownames(summary(model)$coefficients)){
        model.coeffs <- summary(model)$coefficients
        effect.list[i] <- model.coeffs['perturbation', 'Estimate']
        pvalue.list[i] <- model.coeffs['perturbation', 'Pr(>|z|)']
    }
}

# write model outputs to file
model.table <- cbind(
    enhancer.list,
    gene.list,
    effect.list,
    pvalue.list
)
write.csv(
    model.table,
    'data/gasperini/processed/baseline_models_mismatch_gene.csv',
    row.names = FALSE,
    quote = FALSE
)
