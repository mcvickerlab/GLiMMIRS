# Author: Karthik Guruvayurappan

library(stats)
library(rhdf5)
library(MASS)

# create directory for outputs
dir.create('data/experimental/processed/permutation_test_results/')

# read in at-scale enhancer pair analysis results
model.results <- data.frame()

for (i in 1:32) {
    batch.file <- paste0(
        'data/experimental/processed/enhancer_pairs_at_scale_',
        i,
        '.csv'
    )
    batch.results <- read.csv(batch.file)
    model.results <- rbind(model.results, batch.results)
}

# there should be 82,314 model results in total

# filter for results with valid guide efficiency info
model.results <- model.results[complete.cases(model.results), ]

# should have 69,660 total models

# compute FDR-adjusted p-values and filter
model.results$adj.interaction.pvalues <- p.adjust(
    model.results$interaction.pvalues,
    method = 'fdr'
)
significant.results <- model.results[
    model.results$adj.interaction.pvalues < 0.1,

]

# should have 46 significant interactions

# filter for necessary columns in significant interactions
significant.results <- significant.results[
    ,
    c('enhancer.1.list', 'enhancer.2.list', 'gene.list')
]
colnames(significant.results) <- c('enhancer.1', 'enhancer.2', 'gene')


h5.name <- 'data/experimental/processed/gasperini_data.h5'

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

# iterate through significant results and generate permutations
for (i in 1:nrow(significant.results)) {

    # define enhancers and gene
    enhancer.1 <- significant.results[i, 'enhancer.1']
    enhancer.2 <- significant.results[i, 'enhancer.2']
    gene <- significant.results[i, 'gene']

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
    enh.1.no.perturbation <- rep(1, ncol(guide.matrix))
    for (j in 1:nrow(enh.1.efficiencies)) {
        spacer <- enh.1.efficiencies[j, 'spacer']
        efficiency <- enh.1.efficiencies[j, 'Cutting.Efficiency']
        spacer.vector <- guide.matrix[spacer, ]
        spacer.perturbation <- spacer.vector * efficiency
        spacer.no.perturbation <- 1 - spacer.perturbation
        enh.1.no.perturbation <- enh.1.no.perturbation * spacer.no.perturbation
    }
    enh.1.perturbation <- 1 - enh.1.no.perturbation

    enh.2.no.perturbation <- rep(1, ncol(guide.matrix))
    for (j in 1:nrow(enh.2.efficiencies)) {
        spacer <- enh.2.efficiencies[j, 'spacer']
        efficiency <- enh.2.efficiencies[j, 'Cutting.Efficiency']
        spacer.vector <- guide.matrix[spacer, ]
        spacer.perturbation <- spacer.vector * efficiency
        spacer.no.perturbation <- 1 - spacer.perturbation
        enh.2.no.perturbation <- enh.2.no.perturbation * spacer.no.perturbation
    }
    enh.2.perturbation <- 1 - enh.2.no.perturbation

    # get gene counts
    gene.counts <- expr.matrix[gene, ]

    # create data frame and formula for modeling
    model.df <- data.frame(cbind(
        enh.1.perturbation,
        enh.2.perturbation,
        gene.counts,
        covariates
    ))
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

    # initialize empty vectors to store outputs
    permutation.interaction.coeffs <- rep(NA, 100)
    permutation.interaction.pvalues <- rep(NA, 100)

    for (j in 1:100) {

        print(paste0('iteration ', j))

        # shuffle perturbations -> null distribution of interactions
        perturbations.df <- model.df[
            ,
            c('enh.1.perturbation', 'enh.2.perturbation')
        ]
        perturbations.df <- perturbations.df[sample(nrow(perturbations.df)), ]

        model.df$enh.1.perturbation <- perturbations.df$enh.1.perturbation
        model.df$enh.2.perturbation <- perturbations.df$enh.2.perturbation
        
        # refit model with resampled cells
        model <- glm.nb(
            formula = model.formula,
            data = model.df
        )

        # store model outputs
        model.values <- summary(model)$coefficients
        permutation.interaction.coeffs[j] <- model.values[
            'enh.1.perturbation:enh.2.perturbation',
            'Estimate'
        ]
        permutation.interaction.pvalues[j] <- model.values[
            'enh.1.perturbation:enh.2.perturbation',
            'Pr(>|z|)'
        ]
    }

    # write estimates and p-values to outputs
    output.df <- data.frame(cbind(
        permutation.interaction.coeffs,
        permutation.interaction.pvalues
    ))
    write.csv(
        output.df,
        paste0(
            'data/experimental/processed/permutation_test_results/',
            enhancer.1,
            '_',
            enhancer.2,
            '_',
            gene,
            '_permutations.csv'
        ),
        row.names = FALSE,
        quote = FALSE
    )
}
