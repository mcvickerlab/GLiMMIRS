# This program plots a qqplot of the interaction term p-values for the 330 enhancer-enhancer pairs
# determined from the 664 enhancer-gene pairs publihsed by Gasperini et al. (2019). 
#
# Author: Karthik Guruvayurappan

library(stats)
library(ggplot2)
library(RColorBrewer)

# read in at-scale enhancer pair model results across all 32 batches
models <- data.frame()

for (i in 1:32) {
    batch.file.name <- paste0(
        'data/experimental/processed/enhancer_pairs_at_scale_',
        i,
        '.csv'
    )
    batch.models <- read.csv(batch.file.name)
    models <- rbind(models, batch.models)
}

# filter for cases where all coefficients exist
models <- models[complete.cases(models), ]

# pre-process for qqplot
models <- models[order(models$interaction.pvalues), ]
models$adj.interaction.pvalues <- p.adjust(
    models$interaction.pvalues,
    method = 'fdr'
)
models$unif <- seq(1, nrow(models), length.out = nrow(models)) / nrow(models)
models$log.unif <- -log10(models$unif)
models$log.int.pvalue <- -log10(models$interaction.pvalues)
models$is.significant <- models$adj.interaction.pvalues < 0.1

# plot qq-plot
plot <- ggplot(models, aes(
    x = log.unif,
    y = log.int.pvalue,
    color = is.significant
    )) +
    geom_point() +
    theme_classic() + 
    scale_color_manual(
        breaks = c(TRUE, FALSE),
        values = c('red', 'black')
    )

ggsave(
    'out/24_01_24_at_scale_enhancer_pairs_qqplot.png',
    plot,
    device = 'png'
)



# # read in interaction term p-values
# enhancer.enhancer.pvalues <- read.csv('/iblm/netapp/data1/external/Gasperini2019/processed/23_01_12_enhancer_enhancer_at_scale_20_cells_pseudocount_model.csv')
# interaction.pvalues <- enhancer.enhancer.pvalues$interaction.pvalue.list
# interaction.pvalues <- interaction.pvalues[complete.cases(interaction.pvalues)]
# interaction.pvalues <- data.frame(interaction.pvalues[order(interaction.pvalues)])
# colnames(interaction.pvalues) <- c('interaction.pvalue')
# interaction.pvalues$adjusted <- p.adjust(interaction.pvalues$interaction.pvalue, method = 'fdr')
# interaction.pvalues$unif <- seq(1, nrow(interaction.pvalues), length.out = nrow(interaction.pvalues)) / nrow(interaction.pvalues)
# interaction.pvalues$interaction.pvalue <- -log10(interaction.pvalues$interaction.pvalue)
# interaction.pvalues$unif <- -log10(interaction.pvalues$unif)
# interaction.pvalues$color <- '3,808 (insignificant)'
# interaction.pvalues$color[interaction.pvalues$adjusted < 0.1] <- '3,808 (significant)'

# # read in smaller subset p-values
# enhancer.enhancer.pvalues.subset <- read.csv('/iblm/netapp/data1/external/Gasperini2019/processed/23_01_12_enhancer_enhancer_pairs_suppl_table_2_pseudocount_model_enhancer_effects.csv')
# interaction.pvalues.subset <- enhancer.enhancer.pvalues.subset$interaction.pvalue.list
# interaction.pvalues.subset <- interaction.pvalues.subset[complete.cases(interaction.pvalues.subset)]
# interaction.pvalues.subset <- data.frame(interaction.pvalues.subset[order(interaction.pvalues.subset)])
# colnames(interaction.pvalues.subset) <- c('interaction.pvalue')
# interaction.pvalues.subset$adjusted <- p.adjust(interaction.pvalues.subset$interaction.pvalue, method = 'fdr')
# interaction.pvalues.subset$unif <- seq(1, nrow(interaction.pvalues.subset), length.out = nrow(interaction.pvalues.subset)) / nrow(interaction.pvalues.subset)
# interaction.pvalues.subset$interaction.pvalue <- -log10(interaction.pvalues.subset$interaction.pvalue)
# interaction.pvalues.subset$unif <- -log10(interaction.pvalues.subset$unif)
# interaction.pvalues.subset$color <- '330 (insignificant)'
# interaction.pvalues.subset$color[interaction.pvalues.subset$adjusted < 0.1] <- '330 (significant)'

# # read in output results from at-scale enhancer-enhancer analysis
# at.scale.results <- read.csv(
#     '/iblm/netapp/data1/external/Gasperini2019/processed/23_01_12_enhancer_enhancer_at_scale_20_cells_pseudocount_model.csv'
# )
# at.scale.results <- at.scale.results[complete.cases(at.scale.results), ]

# # add FDR adjusted p-values and filter for FDR < 0.1
# at.scale.results$adjusted.interaction.pvalue <- p.adjust(at.scale.results$interaction.pvalue.list, method = 'fdr')
# significant.interactions <- at.scale.results[at.scale.results$adjusted.interaction.pvalue < 0.1, ]
# # filter for necessary columns in significant interactions
# significant.interactions <- significant.interactions[, c('enhancer.1.list', 'enhancer.2.list', 'gene.list', 'interaction.coeff.list', 'interaction.pvalue.list')]
# colnames(significant.interactions) <- c('enhancer.1', 'enhancer.2', 'gene', 'interaction.coeff', 'interaction.pvalue')
# significant.interactions <- significant.interactions[order(significant.interactions$interaction.pvalue), ]

# # read in permutation p-values
# permutation.pvalues <- c()

# for (i in 1:nrow(significant.interactions)) {

#     # read in enhancer and gene names
#     enhancer.1 <- significant.interactions[i, 'enhancer.1']
#     enhancer.2 <- significant.interactions[i, 'enhancer.2']
#     gene <- significant.interactions[i, 'gene']
#     interaction.coefficient <- significant.interactions[i, 'interaction.coeff']
#     coeff.df <- data.frame(interaction.coefficient)
#     colnames(coeff.df) <- c('coeff')

#     # read in file with null distribution of interaction coefficients
#     null.coeffs <- read.csv(paste0('/iblm/netapp/home/karthik/GLiMMIRS/gasperini_data/23_03_31_', enhancer.1, '_', enhancer.2, '_', gene, '_null_interaction_coefficient_estimates.csv'))
#     colnames(null.coeffs) <- c('coeffs')

#     pvalue <- mean(abs(null.coeffs) > interaction.coefficient)
#     permutation.pvalues <- c(permutation.pvalues, pvalue)
# }

# permutation.interactions <- interaction.pvalues[interaction.pvalues$color == '3,808 (significant)', ]
# print(permutation.interactions)
# permutation.interactions$interaction.pvalue <- -log10(permutation.pvalues)
# permutation.interactions$color <- 'Permutation'

# # merge two dataframes together
# interaction.pvalues <- rbind(interaction.pvalues, interaction.pvalues.subset, permutation.interactions)

# qq.plot <- ggplot(interaction.pvalues, aes(x = unif, y = interaction.pvalue, color = color)) +
#     geom_abline(slope = 1, intercept = 0) +
#     geom_point(size = 3) +
#     # geom_abline(slope = 1, intercept = 0) +
#     scale_x_continuous(expand = c(0.02, 0)) +
#     scale_y_continuous(expand = c(0.02, 0)) +
#     xlab(bquote(Expected -log[10](italic(p)))) + 
#     ylab(bquote(Observed -log[10](italic(p)))) +
#     theme_classic() +
#     theme(
#     axis.line = element_line(linewidth = 1),
#     axis.title.x = element_text(size = 20, color = 'black'),
#     axis.title.y = element_text(size = 20, color = 'black'),
#     axis.text = element_text(size = 20, color = 'black'),
#     axis.ticks = element_line(color = 'black', linewidth = 1),
#     axis.ticks.length = unit(2, 'mm'),
#     plot.margin = rep(unit(10, 'mm'), 4),
#     legend.title = element_blank(),
#     legend.position = c(0.28, 0.89),
#     legend.text = element_text(size = 16, color = 'black'),
#     ) +
#     scale_color_manual(values = c('330 (insignificant)' = 'darkgray', '3,808 (insignificant)' = 'black', '3,808 (significant)' = 'red', 'Permutation' = 'blue'))

# ggsave(
#     filename = '/iblm/netapp/home/karthik/GLiMMIRS/plots/23_04_24_interaction_term_qqplot_at_scale.pdf',
#     device = 'pdf',
#     plot = qq.plot
# )

# ggsave(
#     filename = '/iblm/netapp/home/karthik/GLiMMIRS/plots/23_04_24_interaction_term_qqplot_at_scale.png',
#     device = 'png',
#     plot = qq.plot
# )
