# This script plots dotplots of gene expression counts for cells with both
# enhancers perturbed for the 46 significant interactions observed in the
# analysis of 82,314 enhancer pairs.
#
# Author: Karthik Guruvayurappan

library(stats)
library(rhdf5)
library(MASS)
library(ggplot2)

source('src/gasperini/models/analysis_helpers.R')

# create directory for outputs
dir.create('out/interaction_dotplots/')
dir.create('out/prop_zero_plots/')

# get significant interactions from Gasperini analysis
significant_results <- get_significant_results()

# define h5 file name as variable
h5_name <- 'data/gasperini/processed/gasperini_data.h5'

# read in enhancer-enhancer pairs
enhancer_pairs <- h5read(
    h5_name,
    'enhancer_enhancer_at_scale'
)

# read in enhancer-guide table
enhancer_guide <- h5read(
    h5_name,
    'enhancer_guide'
)

# read in cell-guide matrix
guide_matrix <- h5read(
    h5_name,
    'grna/guide_matrix'
)
guide_names <- h5read(
    h5_name,
    'grna/guide_names'
)
rownames(guide_matrix) <- guide_names
barcodes <- h5read(
    h5_name,
    'grna/cell_barcodes'
)
colnames(guide_matrix) <- barcodes

# read in counts matrix
expr_matrix <- h5read(
    h5_name,
    'expr/expr_matrix'
)
genes <- h5read(
    h5_name,
    'expr/gene_names'
)
rownames(expr_matrix) <- genes
barcodes <- h5read(
    h5_name,
    'expr/cell_barcodes'
)
colnames(expr_matrix) <- barcodes

# iterate through significant interactions
for (i in 1:nrow(significant_results)) {

    # get name of enhancers and gene
    enhancer_1 <- significant_results[i, 'enhancer.1.list']
    enhancer_2 <- significant_results[i, 'enhancer.2.list']
    gene <- significant_results[i, 'gene.list']

    # print statement (for progress)
    print(paste0(
        'running ',
        enhancer_1, 
        ' and ', 
        enhancer_2, 
        ' and ',
        gene, '!'
    ))

    # get enhancer spacer sequences for enhancers 1 and 2
    enh_1_spacers <- enhancer_guide[
        enhancer_guide$target.site == enhancer_1,
        'spacer'
    ]
    enh_2_spacers <- enhancer_guide[
        enhancer_guide$target.site == enhancer_2,
        'spacer'
    ]

    # check if guide names are in guide matrix
    enh_1_spacers <- enh_1_spacers[enh_1_spacers %in% guide_names]
    enh_2_spacers <- enh_2_spacers[enh_2_spacers %in% guide_names]

    # compute total perturbation vectors for each enhancer
    enh_1_perturb <- guide_matrix[enh_1_spacers, ]
    enh_2_perturb <- guide_matrix[enh_2_spacers, ]

    if (length(enh_1_spacers) > 1) {
        enh_1_perturb <- colSums(enh_1_perturb)
    }

    if (length(enh_2_spacers) > 1) {
        enh_2_perturb <- colSums(enh_2_perturb)
    }

    # set max perturbation to 1
    enh_1_perturb[enh_1_perturb > 0] <- 1
    enh_2_perturb[enh_2_perturb > 0] <- 1

    # get interaction counts
    interaction_counts <- enh_1_perturb * enh_2_perturb
    interaction_counts <- expr_matrix[gene, interaction_counts > 0]
    interaction_counts <- data.frame(interaction_counts)
    interaction_counts$group <- ''

    # plot interaction dotplot
    plot <- ggplot(interaction_counts, aes(x = group, y = interaction_counts)) + 
        geom_dotplot(binaxis='y', stackdir='center', dotsize = 1.4) +
        theme_classic() +
        xlab('') +
        ylab('Expression') +
        scale_y_continuous(expand = c(0.02, 0)) +
        theme(
          axis.line = element_line(linewidth = 1),
          axis.title.y = element_text(size = 20, color = 'black'),
          axis.text = element_text(size = 20, color = 'black'),
          axis.ticks = element_line(color = 'black', linewidth = 1),
          axis.ticks.length = unit(2, 'mm'),
          plot.margin = rep(unit(10, 'mm'), 4),
          legend.text = element_text(size = 12)
        )

    ggsave(
        paste0(
           'out/interaction_dotplots/',
           enhancer_1,
           '_',
           enhancer_2,
           '_',
           gene,
           '_interaction_dotplot.png' 
        ),
        plot,
        device = 'png'
    )

    # gene.counts <- expr.matrix[gene, ]
    # no.perturbation.counts <- gene.counts[
    #     (enh.1.perturb == 0) & (enh.2.perturb == 0)
    # ]
    # enh.1.perturbation.counts <- gene.counts[
    #     (enh.1.perturb == 1) & (enh.2.perturb == 0)
    # ]
    # enh.2.perturbation.counts <- gene.counts[
    #     (enh.1.perturb == 0) & (enh.2.perturb == 1)
    # ]
    # interaction.counts <- gene.counts[
    #     (enh.1.perturb == 1) & (enh.2.perturb == 1)
    # ]
    # no.perturbation.prop <- sum(no.perturbation.counts == 0) / length(no.perturbation.counts)
    # enh.1.perturbation.prop <- sum(enh.1.perturbation.counts == 0) / length(enh.1.perturbation.counts)
    # enh.2.perturbation.prop <- sum(enh.2.perturbation.counts == 0) / length(enh.2.perturbation.counts)
    # interaction.prop <- sum(interaction.counts == 0) / length(interaction.counts)

    # plot.df <- data.frame(
    #     name = c('No Perturbation', 'E1', 'E2', 'E1 + E2'),
    #     prop = c(no.perturbation.prop, enh.1.perturbation.prop, enh.2.perturbation.prop, interaction.prop)
    # )

    # plot <- ggplot(plot.df, aes(x = name, y = prop)) +
    #     geom_bar(stat = 'identity') +
    #     theme_classic()

    # ggsave(
    #     paste0(
    #        'out/prop_zero_plots/',
    #        enhancer.1,
    #        '_',
    #        enhancer.2,
    #        '_',
    #        gene,
    #        '_prop_zeros_plot.png' 
    #     ),
    #     plot,
    #     device = 'png'
    # )
}
