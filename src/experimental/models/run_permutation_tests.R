# Author: Karthik Guruvayurappan

library(rhdf5)
library(MASS)

# create directory for outputs
dir.create('data/experimental/processed/adaptive_permutation_test_results/')

# read in command line argument (for batch processing)
args = commandArgs(trailingOnly=TRUE)

# define h5 name as a variable
h5.name <- 'data/experimental/processed/gasperini_data.h5'

# read in batch significant interactions
significant.results <- read.csv(
    args[1]
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

# read in guide names (but not matrix)
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

# read in genes (but not counts matrix)
genes <- h5read(
    h5.name,
    'expr/gene_names'
)

# iterate through significant results and generate permutations
for (i in 1:nrow(significant.results)) {

    # define enhancers and gene
    enhancer.1 <- significant.results[i, 'enhancer.1.list']
    enhancer.2 <- significant.results[i, 'enhancer.2.list']
    gene <- significant.results[i, 'gene.list']

    # print statement (for progress tracking)
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

    # implement adaptive permutation scheme
    current.iters <- 100
    not.converged <- TRUE
    permutation.outputs.df <- data.frame()

    while ((current.iters <= 10000) & not.converged) {

        # define empty vectors to hold model outputs
        intercept.effects <- rep(NA, current.iters)
        intercept.pvalues <- rep(NA, current.iters)

        enhancer.1.effects <- rep(NA, current.iters)
        enhancer.1.pvalues <- rep(NA, current.iters)

        enhancer.2.effects <- rep(NA, current.iters)
        enhancer.2.pvalues <- rep(NA, current.iters)

        interaction.effects <- rep(NA, current.iters)
        interaction.pvalues <- rep(NA, current.iters)


        for (j in 1:current.iters) {

            # shuffle perturbations -> null distribution of interactions
            perturbations.df <- model.df[
                ,
                c('enh.1.perturbation', 'enh.2.perturbation')
            ]
            perturbations.df <- perturbations.df[
                sample(nrow(perturbations.df)), 
            ]
            permutation.df <- model.df
            permutation.df$enh.1.perturbation <- perturbations.df$enh.1.perturbation
            permutation.df$enh.2.perturbation <- perturbations.df$enh.2.perturbation

            # refit model with resampled cells
            model <- glm.nb(
                formula = model.formula,
                data = permutation.df
            )

            model.coeffs <- summary(model)$coefficients

            # write out enhancer and interaction coefficients and pvalues
            if ('enh.1.perturbation' %in% rownames(model.coeffs)){
                enhancer.1.effects[j] <- model.coeffs['enh.1.perturbation', 'Estimate']
                enhancer.1.pvalues[j] <- model.coeffs['enh.1.perturbation', 'Pr(>|z|)']
            }
            if ('enh.2.perturbation' %in% rownames(model.coeffs)){
                enhancer.2.effects[j] <- model.coeffs['enh.2.perturbation', 'Estimate']
                enhancer.2.pvalues[j] <- model.coeffs['enh.2.perturbation', 'Pr(>|z|)']
            }
            if ('enh.1.perturbation:enh.2.perturbation' %in% rownames(model.coeffs)){
                interaction.effects[j] <- model.coeffs[
                    'enh.1.perturbation:enh.2.perturbation',
                    'Estimate'
                ]
                interaction.pvalues[j] <- model.coeffs[
                    'enh.1.perturbation:enh.2.perturbation',
                    'Pr(>|z|)'
                ]
            }

            # write intercept information
            intercept.effects[j] <- model.coeffs['(Intercept)', 'Estimate']
            intercept.pvalues[j] <- model.coeffs['(Intercept)', 'Pr(>|z|)']
        }

        true.interaction.coeff <- significant.results[i, 'interaction.effects']
        perm.interaction.coeffs <- interaction.effects[!is.na(interaction.effects)]
        perm.pvalue <- sum(abs(perm.interaction.coeffs) > abs(true.interaction.coeff)) / length(perm.interaction.coeffs)

        # stop early if p-value is insignificant
        if ((perm.pvalue > 0.1) | (current.iters == 10000)) {
            not.converged <- FALSE
            permutation.outputs.df <- data.frame(cbind(
                intercept.effects,
                intercept.pvalues,
                enhancer.1.effects,
                enhancer.1.pvalues,
                enhancer.2.effects,
                enhancer.2.pvalues,
                interaction.effects,
                interaction.pvalues
            ))
        }

        # increase permutations if significant
        current.iters <- current.iters * 10
    }

    # write out permutation estimate and p-values
    write.csv(
        permutation.outputs.df,
        paste0(
            'data/experimental/processed/adaptive_permutation_test_results/',
            enhancer.1,
            '_',
            enhancer.2,
            '_',
            gene,
            '_adaptive_permutations.csv'
        ),
        row.names = FALSE,
        quote = FALSE
    )
}
