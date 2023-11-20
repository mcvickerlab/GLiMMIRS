# This script plots all of the covariate effect sizes for the Gasperini et al.
# data.
#
# Author: Karthik Guruvayurappan

library(ggplot2)

plot.histogram <- function(plot.df, plot.column, var.name) {

    plot <- ggplot(plot.df, aes(x = plot.df[, plot.column])) +
        geom_histogram(color = 'black') +
        scale_x_continuous(expand = c(0, 0)) +
        scale_y_continuous(expand = c(0, 0)) +
        xlab(var.name) + 
        ylab(bquote(Count)) +
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
        filename = paste0('/iblm/netapp/home/karthik/GLiMMIRS/out/23_11_03_', 'at_scale_', plot.column, '.png'),
        device = 'png',
        plot = plot
    )
}

# generate plots for 330 pairs
high.confidence.set <- read.csv('/iblm/netapp/data1/external/Gasperini2019/processed/23_10_31_enhancer_enhancer_pairs_suppl_table_2_pseudocount_model_enhancer_effects.csv')
plot.histogram(high.confidence.set, 'intercept.coeff.list', 'Intercept Effect')
plot.histogram(high.confidence.set, 'enhancer.1.coeff.list', 'Enhancer 1 Effect')
plot.histogram(high.confidence.set, 'enhancer.2.coeff.list', 'Enhancer 2 Effect')
plot.histogram(high.confidence.set, 'interaction.coeff.list', 'Interaction Effect')
plot.histogram(high.confidence.set, 'guide.count.coeff.list', 'gRNA Count Effect')
plot.histogram(high.confidence.set, 'percent.mito.coeff.list', 'Percent Mito Effect')
plot.histogram(high.confidence.set, 's.score.coeff.list', 'S Score Effect')
plot.histogram(high.confidence.set, 'g2m.score.coeff.list', 'G2M Score Effect')

at.scale.set <- read.csv('/iblm/netapp/data1/external/Gasperini2019/processed/23_10_31_enhancer_enhancer_at_scale_20_cells_pseudocount_model.csv')
plot.histogram(at.scale.set, 'intercept.coeff.list', 'Intercept Effect')
plot.histogram(at.scale.set, 'enhancer.1.coeff.list', 'Enhancer 1 Effect')
plot.histogram(at.scale.set, 'enhancer.2.coeff.list', 'Enhancer 2 Effect')
plot.histogram(at.scale.set, 'interaction.coeff.list', 'Interaction Effect')
plot.histogram(at.scale.set, 'guide.count.coeff.list', 'gRNA Count Effect')
plot.histogram(at.scale.set, 'percent.mito.coeff.list', 'Percent Mito Effect')
plot.histogram(at.scale.set, 's.score.coeff.list', 'S Score Effect')
plot.histogram(at.scale.set, 'g2m.score.coeff.list', 'G2M Score Effect')