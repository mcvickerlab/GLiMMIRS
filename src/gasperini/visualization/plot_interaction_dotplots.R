# This script plots dotplots of gene expression counts for cells with both
# enhancers perturbed for the 46 significant interactions observed in the
# analysis of 82,314 enhancer pairs.
#
# Author: Karthik Guruvayurappan

library(stats)
library(rhdf5)
library(MASS)
library(ggplot2)

# create directory for outputs
dir.create('out/interaction_dotplots/')
dir.create('out/prop_zero_plots/')

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

# filter for results with valid guide efficiency info
model.results <- model.results[complete.cases(model.results), ]

# compute FDR-adjusted p-values and filter
model.results$adj.interaction.pvalues <- p.adjust(
    model.results$interaction.pvalues,
    method = 'fdr'
)
significant.results <- model.results[
    model.results$adj.interaction.pvalues < 0.1,

]

# filter for necessary columns in significant interactions
significant.results <- significant.results[
    ,
    c('enhancer.1.list', 'enhancer.2.list', 'gene.list')
]
colnames(significant.results) <- c('enhancer.1', 'enhancer.2', 'gene')

# define h5 file name as variable
h5.name <- 'data/experimental/processed/gasperini_data.h5'

# read in enhancer-enhancer pairs
enhancer.pairs <- h5read(
    h5.name,
    'enhancer_enhancer_at_scale'
)

# read in enhancer-guide table
enhancer.guide <- h5read(
    h5.name,
    'enhancer_guide'
)

# read in cell-guide matrix
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

# iterate through significant interactions
for (i in 1:nrow(significant.results)) {

    # get name of enhancers and gene
    enhancer.1 <- significant.results[i, 'enhancer.1']
    enhancer.2 <- significant.results[i, 'enhancer.2']
    gene <- significant.results[i, 'gene']

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

    # check if guide names are in guide matrix
    enh.1.spacers <- enh.1.spacers[enh.1.spacers %in% guide.names]
    enh.2.spacers <- enh.2.spacers[enh.2.spacers %in% guide.names]

    # compute total perturbation vectors for each enhancer
    enh.1.perturb <- guide.matrix[enh.1.spacers, ]
    enh.2.perturb <- guide.matrix[enh.2.spacers, ]

    if (length(enh.1.spacers) > 1) {
        enh.1.perturb <- colSums(enh.1.perturb)
    }

    if (length(enh.2.spacers) > 1) {
        enh.2.perturb <- colSums(enh.2.perturb)
    }

    # set max perturbation to 1
    enh.1.perturb[enh.1.perturb > 0] <- 1
    enh.2.perturb[enh.2.perturb > 0] <- 1

    # get interaction counts
    interaction.counts <- enh.1.perturb * enh.2.perturb
    interaction.counts <- expr.matrix[gene, interaction.counts > 0]
    interaction.counts <- data.frame(interaction.counts)
    interaction.counts$group <- 'group'

    # plot interaction dotplot
    plot <- ggplot(interaction.counts, aes(x = group, y = interaction.counts)) + 
        geom_dotplot(binaxis='y', stackdir='center', dotsize = 1.6) +
        theme_classic()

    ggsave(
        paste0(
           'out/interaction_dotplots/',
           enhancer.1,
           '_',
           enhancer.2,
           '_',
           gene,
           '_interaction_dotplot.png' 
        ),
        plot,
        device = 'png'
    )

    gene.counts <- expr.matrix[gene, ]
    no.perturbation.counts <- gene.counts[
        (enh.1.perturb == 0) & (enh.2.perturb == 0)
    ]
    enh.1.perturbation.counts <- gene.counts[
        (enh.1.perturb == 1) & (enh.2.perturb == 0)
    ]
    enh.2.perturbation.counts <- gene.counts[
        (enh.1.perturb == 0) & (enh.2.perturb == 1)
    ]
    interaction.counts <- gene.counts[
        (enh.1.perturb == 1) & (enh.2.perturb == 1)
    ]
    no.perturbation.prop <- sum(no.perturbation.counts == 0) / length(no.perturbation.counts)
    enh.1.perturbation.prop <- sum(enh.1.perturbation.counts == 0) / length(enh.1.perturbation.counts)
    enh.2.perturbation.prop <- sum(enh.2.perturbation.counts == 0) / length(enh.2.perturbation.counts)
    interaction.prop <- sum(interaction.counts == 0) / length(interaction.counts)

    plot.df <- data.frame(
        name = c('No Perturbation', 'E1', 'E2', 'E1 + E2'),
        prop = c(no.perturbation.prop, enh.1.perturbation.prop, enh.2.perturbation.prop, interaction.prop)
    )

    plot <- ggplot(plot.df, aes(x = name, y = prop)) +
        geom_bar(stat = 'identity') +
        theme_classic()

    ggsave(
        paste0(
           'out/prop_zero_plots/',
           enhancer.1,
           '_',
           enhancer.2,
           '_',
           gene,
           '_prop_zeros_plot.png' 
        ),
        plot,
        device = 'png'
    )
}
