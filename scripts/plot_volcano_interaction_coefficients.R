# This program plots a volcano plots of the interaction coefficients from the at-scale interaction
# analysis and their associated p-values.
#
# Author: Karthik Guruvayurappan

library(stats)
library(BoutrosLab.plotting.general)

# read in interaction term p-values
enhancer.enhancer.pvalues <- read.csv('/iblm/netapp/data1/external/Gasperini2019/processed/enhancer_enhancer_at_scale_20_cells_pseudocount_model.csv')

# convert pvalues to -log10 scale
interaction.pvalues <- enhancer.enhancer.pvalues[, c('interaction.coeff.list', 'interaction.pvalue.list')]
colnames(interaction.pvalues) <- c('coeff', 'pvalue')
interaction.pvalues$scaled.pvalue <- -1 * log10(interaction.pvalues$pvalue)
interaction.pvalues$adjusted.pvalue <- p.adjust(interaction.pvalues$pvalue, method = 'fdr')
interaction.pvalues$dot.color <- 'black'
interaction.pvalues$dot.color[interaction.pvalues$adjusted.pvalue < 0.1] <- 'red'

create.scatterplot(
    formula = scaled.pvalue ~ coeff,
    data = interaction.pvalues,
    filename = '/iblm/netapp/home/karthik/crisprQTL/plots/at_scale_interaction_volcano_plot.tiff',
    resolution = 300,
    xlimits = c(-15, 15),
    col = interaction.pvalues$dot.color,
    xlab.label = 'Interaction',
    ylab.label = '-log10(pvalue)'
)


