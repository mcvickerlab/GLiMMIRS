# This script runs the GLiMMIRS model on the at-scale 82,314 enhancer-enhancer
# pairs derived from the Gasperini et al. data.
#
# Author: Karthik Guruvayurappan

library(rhdf5)
library(MASS)

# read in command line argument (for batch processing)
args = commandArgs(trailingOnly=TRUE)

# # define h5 file name as a variable
h5.name <- 'data/experimental/processed/gasperini_data.h5'

# read in enhancer-enhancer pairs
enhancer.pairs <- read.csv(
    args[1]
)

# read in enhancer-guide table
enhancer.guide <- h5read(
    h5.name,
    'enhancer_guide'
)

# # read in guide efficiency information
guide.info <- h5read(
    h5.name,
    'grna/guide_info'
)

# # read in guide matrix info but not matrix (for memory)
# guide.matrix <- h5read(
#     h5.name,
#     'grna/guide_matrix'
# )
guide.names <- h5read(
    h5.name,
    'grna/guide_names'
)

# filter enhancer-guide table to only include guides in matrix
enhancer.guide <- enhancer.guide[enhancer.guide$spacer %in% guide.names, ]

# read in cell covariates
covariates <- h5read(
    h5.name,
    'expr/cell_covariates'
)

# # read in genes but not counts matrix (for memory)
# expr.matrix <- h5read(
#     h5.name,
#     'expr/expr_matrix'
# )
genes <- h5read(
    h5.name,
    'expr/gene_names'
)
# rownames(expr.matrix) <- genes
# barcodes <- h5read(
#     h5.name,
#     'expr/cell_barcodes'
# )
# colnames(expr.matrix) <- barcodes

# # create vectors to hold model outputs
enhancer.1.list <- rep(NA, nrow(enhancer.pairs))
enhancer.2.list <- rep(NA, nrow(enhancer.pairs))
gene.list <- rep(NA, nrow(enhancer.pairs))

intercept.effects <- rep(NA, nrow(enhancer.pairs))
intercept.pvalues <- rep(NA, nrow(enhancer.pairs))

enhancer.1.effects <- rep(NA, nrow(enhancer.pairs))
enhancer.1.pvalues <- rep(NA, nrow(enhancer.pairs))

enhancer.2.effects <- rep(NA, nrow(enhancer.pairs))
enhancer.2.pvalues <- rep(NA, nrow(enhancer.pairs))

interaction.effects <- rep(NA, nrow(enhancer.pairs))
interaction.pvalues <- rep(NA, nrow(enhancer.pairs))

prep.batch.effects <- rep(NA, nrow(enhancer.pairs))
prep.batch.pvalues <- rep(NA, nrow(enhancer.pairs))

guide.count.effects <- rep(NA, nrow(enhancer.pairs))
guide.count.pvalues <- rep(NA, nrow(enhancer.pairs))

percent.mito.effects <- rep(NA, nrow(enhancer.pairs))
percent.mito.pvalues <- rep(NA, nrow(enhancer.pairs))

s.score.effects <- rep(NA, nrow(enhancer.pairs))
s.score.pvalues <- rep(NA, nrow(enhancer.pairs))

g2m.score.effects <- rep(NA, nrow(enhancer.pairs))
g2m.score.pvalues <- rep(NA, nrow(enhancer.pairs))


for (i in 1:nrow(enhancer.pairs)) {

    # get name of enhancers and gene
    enhancer.1 <- enhancer.pairs[i, 'enhancer_1']
    enhancer.2 <- enhancer.pairs[i, 'enhancer_2']
    gene <- enhancer.pairs[i, 'gene']

    # print statement (for progress)
    print(paste0(
        'running ',
        enhancer.1, 
        ' and ', 
        enhancer.2, 
        ' and ',
        gene, '!'
    ))

    # get enhancer spacer sequences for enhancers 1 and 2
    enh.1.spacers <- enhancer.guide[
        enhancer.guide$target.site == enhancer.1,
        'spacer'
    ]
    enh.2.spacers <- enhancer.guide[
        enhancer.guide$target.site == enhancer.2,
        'spacer'
    ]

    # get spacer guide efficiencies and set NAs to 0
    enh.1.efficiencies <- guide.info[guide.info$spacer %in% enh.1.spacers, ]
    enh.1.efficiencies <- enh.1.efficiencies[
        , 
        c('spacer', 'Cutting.Efficiency')
    ]
    enh.1.efficiencies[is.na(enh.1.efficiencies)] <- 0

    enh.2.efficiencies <- guide.info[guide.info$spacer %in% enh.2.spacers, ]
    enh.2.efficiencies <- enh.2.efficiencies[
        , 
        c('spacer', 'Cutting.Efficiency')
    ]
    enh.2.efficiencies[is.na(enh.2.efficiencies)] <- 0

    # compute guide perturbation vector by computing 1 - (no perturbation prob)
    enh.1.no.perturbation <- rep(1, nrow(covariates))
    for (j in 1:nrow(enh.1.efficiencies)) {
        spacer <- enh.1.efficiencies[j, 'spacer']
        efficiency <- enh.1.efficiencies[j, 'Cutting.Efficiency']
        spacer.index <- match(spacer, guide.names)
        spacer.vector <- h5read(
            h5.name,
            'grna/guide_matrix',
            index = list(spacer.index, NULL)
        )
        spacer.perturbation <- spacer.vector * efficiency
        spacer.no.perturbation <- 1 - spacer.perturbation
        enh.1.no.perturbation <- enh.1.no.perturbation * spacer.no.perturbation
    }
    enh.1.perturbation <- 1 - enh.1.no.perturbation

    enh.2.no.perturbation <- rep(1, nrow(covariates))
    for (j in 1:nrow(enh.2.efficiencies)) {
        spacer <- enh.2.efficiencies[j, 'spacer']
        efficiency <- enh.2.efficiencies[j, 'Cutting.Efficiency']
        spacer.index <- match(spacer, guide.names)
        spacer.vector <- h5read(
            h5.name,
            'grna/guide_matrix',
            index = list(spacer.index, NULL)
        )
        spacer.perturbation <- spacer.vector * efficiency
        spacer.no.perturbation <- 1 - spacer.perturbation
        enh.2.no.perturbation <- enh.2.no.perturbation * spacer.no.perturbation
    }
    enh.2.perturbation <- 1 - enh.2.no.perturbation

    # get gene counts
    gene.index <- match(gene, genes)
    gene.counts <- h5read(
        h5.name,
        'expr/expr_matrix',
        index = list(gene.index, NULL)
    )

    # add pseudocount to gene counts
    pseudocount <- 0.01
    gene.counts <- gene.counts + pseudocount

    # fit negative binomial GLM
    enh.1.perturbation <- t(enh.1.perturbation)[,1]
    enh.2.perturbation <- t(enh.2.perturbation)[,1]
    gene.counts <- t(gene.counts)[,1]
    covariates$scaling.factor <- as.vector(covariates$scaling.factor)
    model.df <- cbind(
        enh.1.perturbation,
        enh.2.perturbation,
        gene.counts,
        covariates
    )
    model.formula <- as.formula(paste0(
        'gene.counts ~ ',
        'enh.1.perturbation * ',
        'enh.2.perturbation + ',
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
    enhancer.1.list[i] <- enhancer.1
    enhancer.2.list[i] <- enhancer.2
    gene.list[i] <- gene

    model.coeffs <- summary(model)$coefficients

    # write out enhancer and interaction coefficients and pvalues
    if ('enh.1.perturbation' %in% rownames(model.coeffs)){
        enhancer.1.effects[i] <- model.coeffs['enh.1.perturbation', 'Estimate']
        enhancer.1.pvalues[i] <- model.coeffs['enh.1.perturbation', 'Pr(>|z|)']
    }
    if ('enh.2.perturbation' %in% rownames(model.coeffs)){
        enhancer.2.effects[i] <- model.coeffs['enh.2.perturbation', 'Estimate']
        enhancer.2.pvalues[i] <- model.coeffs['enh.2.perturbation', 'Pr(>|z|)']
    }
    if ('enh.1.perturbation:enh.2.perturbation' %in% rownames(model.coeffs)){
        interaction.effects[i] <- model.coeffs[
            'enh.1.perturbation:enh.2.perturbation',
            'Estimate'
        ]
        interaction.pvalues[i] <- model.coeffs[
            'enh.1.perturbation:enh.2.perturbation',
            'Pr(>|z|)'
        ]
    }

    # write intercept information
    intercept.effects[i] <- model.coeffs['(Intercept)', 'Estimate']
    intercept.pvalues[i] <- model.coeffs['(Intercept)', 'Pr(>|z|)']

    # store covariate information
    prep.batch.effects[i] <- model.coeffs[
        'prep_batchprep_batch_2',
        'Estimate'
    ]
    prep.batch.pvalues[i] <- model.coeffs[
        'prep_batchprep_batch_2',
        'Pr(>|z|)'
    ]

    guide.count.effects[i] <- model.coeffs['guide_count', 'Estimate']
    guide.count.pvalues[i] <- model.coeffs['guide_count', 'Pr(>|z|)']

    percent.mito.effects[i] <- model.coeffs['percent.mito', 'Estimate']
    percent.mito.pvalues[i] <- model.coeffs['percent.mito', 'Pr(>|z|)']

    s.score.effects[i] <- model.coeffs['s.score', 'Estimate']
    s.score.pvalues[i] <- model.coeffs['s.score', 'Pr(>|z|)']

    g2m.score.effects[i] <- model.coeffs['g2m.score', 'Estimate']
    g2m.score.pvalues[i] <- model.coeffs['g2m.score', 'Pr(>|z|)']
}

# write model outputs to file
model.table <- cbind(
    enhancer.1.list,
    enhancer.2.list,
    gene.list,
    intercept.effects,
    intercept.pvalues,
    enhancer.1.effects,
    enhancer.1.pvalues,
    enhancer.2.effects,
    enhancer.2.pvalues,
    interaction.effects,
    interaction.pvalues,
    prep.batch.effects,
    prep.batch.pvalues,
    guide.count.effects,
    guide.count.pvalues,
    percent.mito.effects,
    percent.mito.pvalues,
    s.score.effects,
    s.score.pvalues,
    g2m.score.effects,
    g2m.score.pvalues
)

batch.number <- strsplit(args[1], '_')[[1]][3]
batch.number <- strsplit(batch.number, '.', fixed = TRUE)[[1]][1]
write.csv(
    model.table,
    paste0(
        'data/experimental/processed/enhancer_pairs_at_scale_',
        batch.number,
        '.csv'
    ),
    row.names = FALSE,
    quote = FALSE
)
