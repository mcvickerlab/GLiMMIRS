# This program plots a volcano plots of the interaction coefficients from the 
# at-scale interaction analysis and their associated p-values.
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

# perform multiple testing correction on the interaction term p-values
models$adj_interaction_pvalues <- p.adjust(
  models$interaction.pvalues,
  method = 'fdr'
)

# add a variable which marks significant interactions
models$is_significant <- models$adj_interaction_pvalues < 0.1

# plot volcano plot
plot <- ggplot(models, aes(
  x = interaction.effects, 
  y = -log10(interaction.pvalues),
  color = is_significant
  )) +
  geom_point(size = 3) +
  theme_classic() +
  scale_x_continuous(expand = c(0.02, 0)) +
  scale_y_continuous(expand = c(0.02, 0)) +
  xlab('Interaction Coefficient') +
  ylab(bquote(Observed -log[10](italic(p)))) +
  theme(
    axis.line = element_line(linewidth = 1),
    axis.title.x = element_text(size = 24, color = 'black'),
    axis.title.y = element_text(size = 24, color = 'black'),
    axis.text = element_text(size = 24, color = 'black'),
    axis.ticks = element_line(color = 'black', linewidth = 1),
    axis.ticks.length = unit(2, 'mm'),
    plot.margin = rep(unit(10, 'mm'), 4),
    legend.position = 'none'
    # legend.text = element_text(size = 12)
  ) +
  labs(color = NULL) +
  scale_color_manual(
    breaks = c(TRUE, FALSE),
    values = c('red', 'darkgray')
    # labels = c('Significant (FDR < 0.1)', 'Insignificant')
  )

# save plot to output file
ggsave(
  filename = 'out/volcano_interaction_terms_at_scale.pdf',
  plot = plot,
  device = 'pdf',
  width = 7,
  height = 7,
  unit = 'in'
)
