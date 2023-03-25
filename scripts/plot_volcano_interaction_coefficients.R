# This program plots a volcano plots of the interaction coefficients from the at-scale interaction
# analysis and their associated p-values.
#
# Author: Karthik Guruvayurappan

library(stats)
library(ggplot2)
library(RColorBrewer)
# library(BoutrosLab.plotting.general)

# read in interaction term p-values
enhancer.enhancer.pvalues <- read.csv('/iblm/netapp/data1/external/Gasperini2019/processed/23_01_12_enhancer_enhancer_at_scale_20_cells_pseudocount_model.csv')

# convert pvalues to -log10 scale
interaction.pvalues <- enhancer.enhancer.pvalues[, c('interaction.coeff.list', 'interaction.pvalue.list')]
colnames(interaction.pvalues) <- c('coeff', 'pvalue')
interaction.pvalues$scaled.pvalue <- -1 * log10(interaction.pvalues$pvalue)
interaction.pvalues$adjusted.pvalue <- p.adjust(interaction.pvalues$pvalue, method = 'fdr')
interaction.pvalues$dot.color <- 'black'
interaction.pvalues$dot.color[interaction.pvalues$adjusted.pvalue < 0.1] <- 'red'
interaction.pvalues <- interaction.pvalues[complete.cases(interaction.pvalues$coeff), ]

qq.plot <- ggplot(interaction.pvalues, aes(x = coeff, y = scaled.pvalue)) +
    geom_point(aes(color = dot.color)) +
    scale_x_continuous(expand = c(0.02, 0)) +
    scale_y_continuous(expand = c(0.02, 0)) +
    xlab(bquote(Coefficient)) + 
    ylab(bquote(-log[10](italic(p)))) +
    theme_classic() +
    theme(
    axis.line = element_line(linewidth = 1),
    axis.title.x = element_text(size = 20, color = 'black'),
    axis.title.y = element_text(size = 20, color = 'black'),
    axis.text = element_text(size = 20, color = 'black'),
    axis.ticks = element_line(color = 'black', linewidth = 1),
    axis.ticks.length = unit(2, 'mm'),
    plot.margin = rep(unit(10, 'mm'), 4),
    legend.position = "none"
    ) +
    scale_color_manual(values = c('black', 'red'))
    
ggsave(
    filename = '/iblm/netapp/home/karthik/GLiMMIRS/plots/23_03_25_at_scale_interaction_volcano_plot.pdf',
    device = 'pdf',
    plot = qq.plot
)

# create.scatterplot(
#     formula = scaled.pvalue ~ coeff,
#     data = interaction.pvalues,
#     filename = '/iblm/netapp/home/karthik/crisprQTL/plots/23_02_21_at_scale_interaction_volcano_plot.pdf',
#     resolution = 300,
#     xlimits = c(-15, 15),
#     col = interaction.pvalues$dot.color,
#     xlab.label = 'Interaction',
#     ylab.label = '-log10(pvalue)'
# )


