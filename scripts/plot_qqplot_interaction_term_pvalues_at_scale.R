# This program plots a qqplot of the interaction term p-values for the 330 enhancer-enhancer pairs
# determined from the 664 enhancer-gene pairs publihsed by Gasperini et al. (2019). 
#
# Author: Karthik Guruvayurappan

library(stats)
library(BoutrosLab.plotting.general)

# read in interaction term p-values
enhancer.enhancer.pvalues <- read.csv('/iblm/netapp/data1/external/Gasperini2019/processed/23_01_12_enhancer_enhancer_at_scale_20_cells_pseudocount_model.csv')
interaction.pvalues <- enhancer.enhancer.pvalues$interaction.pvalue.list
interaction.pvalues <- interaction.pvalues[complete.cases(interaction.pvalues)]
interaction.pvalues <- data.frame(interaction.pvalues[order(interaction.pvalues)])
print(nrow(interaction.pvalues))
colnames(interaction.pvalues) <- c('interaction.pvalue')
interaction.pvalues$adjusted <- p.adjust(interaction.pvalues$interaction.pvalue, method = 'fdr')
print(sum(interaction.pvalues$adjusted < 0.1))
interaction.pvalues$unif <- seq(1, nrow(interaction.pvalues), length.out = nrow(interaction.pvalues)) / nrow(interaction.pvalues)
interaction.pvalues$interaction.pvalue <- -log10(interaction.pvalues$interaction.pvalue)
interaction.pvalues$unif <- -log10(interaction.pvalues$unif)
interaction.pvalues$color <- 'black'
interaction.pvalues$color[interaction.pvalues$adjusted < 0.1] <- 'red'

create.scatterplot(
    formula = interaction.pvalue ~ unif,
    data = interaction.pvalues,
    col = interaction.pvalues$color,
    filename = '/iblm/netapp/home/karthik/crisprQTL/plots/23_02_21_interaction_term_at_scale_qqplot.pdf',
    resolution = 300,
    add.xyline = TRUE,
    xlab.label = '-log10(expected pvalue)',
    ylab.label = '-log10(observed pvalue)'
)

sorted.df <- enhancer.enhancer.pvalues[order(enhancer.enhancer.pvalues$interaction.pvalue.list), ]
head(sorted.df)
