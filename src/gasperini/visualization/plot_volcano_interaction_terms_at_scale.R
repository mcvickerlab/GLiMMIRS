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
  geom_point() +
  theme_classic() +
  scale_x_continuous(expand = c(0.02, 0)) +
  scale_y_continuous(expand = c(0.02, 0)) +
  xlab('Interaction Coefficient') +
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

# save plot to output file
ggsave(
  filename = 'out/volcano_interaction_terms_at_scale.png',
  plot = plot,
  device = 'png',
  width = 8,
  height = 7,
  unit = 'in'
)




# convert pvalues to -log10 scale
interaction.pvalues$adjusted.pvalue <- p.adjust(interaction.pvalues$pvalue, method = 'fdr')
interaction.pvalues$dot.color <- 'black'
interaction.pvalues$dot.color[interaction.pvalues$adjusted.pvalue < 0.1] <- 'red'
interaction.pvalues <- interaction.pvalues[complete.cases(interaction.pvalues$coeff), ]

volcano.plot <- ggplot(interaction.pvalues, aes(x = coeff, y = scaled.pvalue)) +
    geom_point(aes(color = dot.color), size = 3) +
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
    filename = '/iblm/netapp/home/karthik/GLiMMIRS/plots/23_04_20_at_scale_interaction_volcano_plot.pdf',
    device = 'pdf',
    plot = volcano.plot
)

ggsave(
    filename = '/iblm/netapp/home/karthik/GLiMMIRS/plots/23_04_20_at_scale_interaction_volcano_plot.png',
    device = 'png',
    plot = volcano.plot
)
