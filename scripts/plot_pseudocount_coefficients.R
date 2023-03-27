# This program plots the interaction coefficients obtained when using a pseudocount model against
# the coefficients obtained when using a model without a pseudocount. 
# This program was written by Karthik Guruvayurappan.

library(ggplot2)
library(RColorBrewer)

no.pseudocount.pvalues <- read.csv('/iblm/netapp/data1/external/Gasperini2019/processed/23_03_27_enhancer_enhancer_pairs_suppl_table_2_no_pseudocount_model_enhancer_effects.csv')

pseudocount.pvalues <- read.csv('/iblm/netapp/data1/external/Gasperini2019/processed/23_03_27_enhancer_enhancer_pairs_suppl_table_2_pseudocount_model_enhancer_effects.csv')

no.pseudocount.coefficients <- no.pseudocount.pvalues$interaction.coeff.list
pseudocount.coefficients <- pseudocount.pvalues$interaction.coeff.list
coefficient.df <- data.frame(cbind(pseudocount.coefficients, no.pseudocount.coefficients))

qq.plot <- ggplot(coefficient.df, aes(x = no.pseudocount.coefficients, y = pseudocount.coefficients)) + 
    geom_point(color = 'black') +
    geom_abline(slope = 1, intercept = 0) +
    scale_x_continuous(expand = c(0.02, 0)) +
    scale_y_continuous(expand = c(0.02, 0)) +
    xlab(bquote('No Pseudocount Interaction Coefficient')) + 
    ylab(bquote('Pseudocount Interaction Coefficient')) +
    theme_classic() +
    theme(
        axis.line = element_line(linewidth = 1),
        axis.title.x = element_text(size = 20, color = 'black'),
        axis.title.y = element_text(size = 20, color = 'black'),
        axis.text = element_text(size = 20, color = 'black'),
        axis.ticks = element_line(color = 'black', linewidth = 1),
        axis.ticks.length = unit(2, 'mm'),
        plot.margin = rep(unit(10, 'mm'), 4),
    )

ggsave(
    filename = '/iblm/netapp/home/karthik/GLiMMIRS/plots/23_03_27_pseudocount_no_pseudocount_interaction_coefficients.pdf',
    device = 'pdf',
    plot = qq.plot
)
