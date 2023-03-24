# This program plots a qqplot of the interaction term p-values for the at-scale enhancer-enhancer
# publihsed by Gasperini et al. (2019). 
#
# Author: Karthik Guruvayurappan

library(stats)
library(ggplot2)
library(RColorBrewer)

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

qq.plot <- ggplot(interaction.pvalues, aes(x = unif, y = interaction.pvalue)) +
    geom_point() +
    geom_abline(slope = 1, intercept = 0) +
    scale_x_continuous(expand = c(0.02, 0)) +
    scale_y_continuous(expand = c(0.02, 0)) +
    xlab(bquote(Expected -log[10](italic(p)))) + 
    ylab(bquote(Observed -log[10](italic(p)))) +
    theme_classic() +
    theme(
    axis.line = element_line(linewidth = 1),
    axis.title.x = element_text(size = 20, color = 'black'),
    axis.title.y = element_text(size = 20, color = 'black'),
    axis.text = element_text(size = 20, color = 'black'),
    axis.ticks = element_line(color = 'black', linewidth = 1),
    axis.ticks.length = unit(2, 'mm'),
    plot.margin = rep(unit(10, 'mm'), 4)
    )

ggsave(
    filename = '/iblm/netapp/home/karthik/GLiMMIRS/plots/23_03_24_interaction_term_qqplot.tiff',
    device = 'tiff',
    plot = qq.plot
)
