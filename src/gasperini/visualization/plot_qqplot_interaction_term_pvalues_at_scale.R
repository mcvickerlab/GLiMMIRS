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
        'data/gasperini/processed/enhancer_pairs_at_scale_',
        i,
        '.csv'
    )
    batch.models <- read.csv(batch.file.name)
    models <- rbind(models, batch.models)
}

# filter for cases where all coefficients exist
models <- models[complete.cases(models), ]

# print out how many models remain after filtering
print(dim(models))

# sort results by p-value
models <- models[order(models$interaction.pvalues), ]

# add a uniform distribution (for plotting)
models$unif <- 1:nrow(models) / nrow(models)

# perform multiple testing correction on interaction p-values
models$adj_interaction_pvalues <- p.adjust(
    models$interaction.pvalues,
    method = 'fdr'
)
print(sum(models$adj_interaction_pvalues < 0.1))

# add whether interaction is significant after correction to models
models$is_significant <- models$adj_interaction_pvalues < 0.1

# plot qq-plot
plot <- ggplot(models, aes(
    x = -log10(unif),
    y = -log10(interaction.pvalues),
    color = is_significant
    )) +
    geom_point(size = 2) +
    geom_abline(slope = 1, intercept = 0) +
    theme_classic() + 
    scale_x_continuous(expand = c(0.02, 0)) +
    scale_y_continuous(expand = c(0.02, 0)) +
    xlab(bquote(Expected -log[10](italic(p)))) + 
    ylab(bquote(Observed -log[10](italic(p)))) +
    theme(
      axis.line = element_line(linewidth = 1),
      axis.title.x = element_text(size = 20, color = 'black'),
      axis.title.y = element_text(size = 20, color = 'black'),
      axis.text = element_text(size = 20, color = 'black'),
      axis.ticks = element_line(color = 'black', linewidth = 1),
      axis.ticks.length = unit(2, 'mm'),
      plot.margin = rep(unit(10, 'mm'), 4),
      legend.text = element_text(size = 12)
    ) +
    labs(color = NULL) +
    scale_color_manual(
        breaks = c(TRUE, FALSE),
        values = c('red', 'black'),
        labels = c('Significant (FDR < 0.1)', 'Insignificant')
    )

ggsave(
    'out/qqplot_interaction_pvalues_at_scale.png',
    plot,
    device = 'png',
    width = 8,
    height = 7,
    units = 'in'
)
