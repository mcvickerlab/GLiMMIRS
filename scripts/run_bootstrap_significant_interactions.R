# From the at-scale enhancer-enhancer interaction analysis, there were 6
# interaction terms which were significant. However, it appears that these
# results were dominated by a high count from a single cell. This script
# generates 95% confidence intervals using a boostrap analysis as opposed to
# standard error estimates from the lm() function.
#
# Author: Karthik Guruvayurappan

library(stats)
library(rhdf5)
library(MASS)

# read in output results from at-scale enhancer-enhancer analysis
at.scale.results <- read.csv(
    '/iblm/netapp/data1/external/Gasperini2019/processed/enhancer_enhancer_at_scale_20_cells_pseudocount_model.csv'
)
at.scale.results <- at.scale.results[complete.cases(at.scale.results), ]

# add FDR adjusted p-values and filter for FDR < 0.1
at.scale.results$adjusted.interaction.pvalue <- p.adjust(at.scale.results$interaction.pvalue.list, method = 'fdr')
significant.interactions <- at.scale.results[at.scale.results$adjusted.interaction.pvalue < 0.1, ]

# filter for necessary columns in significant interactions
significant.interactions <- significant.interactions[, c('enhancer.1.list', 'enhancer.2.list', 'gene.list')]
colnames(significant.interactions) <- c('enhancer.1', 'enhancer.2', 'gene')

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

for (i in 1:nrow(significant.interactions)) {

    # define enhancers and gene
    enhancer.1 <- significant.interactions[i, 'enhancer.1']
    enhancer.2 <- significant.interactions[i, 'enhancer.2']
    gene <- significant.interactions[i, 'gene']

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
    pseudocount <- 0.01
    gene.counts <- counts.matrix[gene, ] + pseudocount

    # create dataframe for modeling
    model.df <- cbind(covariates, enhancer.1.indicator.vector, enhancer.2.indicator.vector, gene.counts)

    interaction.coefficient.estimates <- rep(NA, 100)

    for (j in 1:100) {

        print(paste0('iteration ', j))

        # resample cells with replacement
        bootstrap.df <- model.df[sample(1:nrow(model.df), size = nrow(model.df), replace = TRUE), ]
        
        # refit model with resampled cells
        bootstrap.model <- glm.nb(
            formula = gene.counts ~ enhancer.1.indicator.vector * enhancer.2.indicator.vector + prep_batch + guide_count + percent.mito + s.score + g2m.score + offset(log(scaling.factors)),
            data = bootstrap.df
        )

        # record coefficient estimate
        interaction.coefficient.estimates[j] <- summary(bootstrap.model)$coefficients['enhancer.1.indicator.vector:enhancer.2.indicator.vector', 'Estimate']

    }

    # write coefficient estimates to output file
    write.csv(interaction.coefficient.estimates, paste0('/iblm/netapp/data1/external/Gasperini2019/processed/', '22_12_21_', enhancer.1, '_', enhancer.2, '_', gene, '_bootstrap_coefficient_estimates.csv'), row.names = FALSE)
}
