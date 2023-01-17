# This program plots a qqplot of the interaction term p-values for the at-scale enhancer-enhancer
# publihsed by Gasperini et al. (2019). 
#
# Author: Karthik Guruvayurappan

library(stats)
library(BoutrosLab.plotting.general)

# read in interaction term p-values
enhancer.enhancer.pvalues <- read.csv('/iblm/netapp/data1/external/Gasperini2019/processed/23_01_12_enhancer_enhancer_pairs_suppl_table_2_pseudocount_model_enhancer_effects.csv')
interaction.pvalues <- enhancer.enhancer.pvalues$interaction.pvalue.list
interaction.pvalues <- interaction.pvalues[complete.cases(interaction.pvalues)]
interaction.pvalues <- data.frame(interaction.pvalues[order(interaction.pvalues)])
colnames(interaction.pvalues) <- c('interaction.pvalue')
interaction.pvalues$adjusted.pvalue <- p.adjust(interaction.pvalues$interaction.pvalue, method = 'fdr')
interaction.pvalues$unif <- seq(1, nrow(interaction.pvalues), length.out = nrow(interaction.pvalues)) / nrow(interaction.pvalues)
interaction.pvalues$interaction.pvalue <- -log10(interaction.pvalues$interaction.pvalue)
interaction.pvalues$unif <- -log10(interaction.pvalues$unif)
interaction.pvalues$color <- 'black'
interaction.pvalues$color[interaction.pvalues$adjusted.pvalue < 0.1] <- 'red'

create.scatterplot(
    formula = interaction.pvalue ~ unif,
    data = interaction.pvalues,
    col = interaction.pvalues$color,
    filename = '/iblm/netapp/home/karthik/crisprQTL/plots/23_01_13_interaction_term_qqplot.tiff',
    resolution = 300,
    add.xyline = TRUE,
    xlab.label = '-log10(expected pvalue)',
    ylab.label = '-log10(observed pvalue)',
)

sorted.df <- enhancer.enhancer.pvalues[order(enhancer.enhancer.pvalues$interaction.pvalue.list), ]
head(sorted.df)


