# This script plots the distribution of interaction effect sizes and the
# distribution of enhancer effect sizes for the 664 enhancer-gene pairs.
#
# Author: Karthik Guruvayurappan

source('src/gasperini/models/analysis_helpers.R')

library(stats)
library(ggplot2)
library(RColorBrewer)

# read in significant results
significant_results <- get_significant_results()

# isolate interaction effects
interaction_effects <- significant_results$interaction.effects

# read in results from Gasperini et al. re-analysis
gasperini_results <- read.csv(
  'data/gasperini/processed/baseline_models.csv'
)

# filter for valid effects
gasperini_results <- gasperini_results[complete.cases(gasperini_results), ]

# isolate enhancer effects
enhancer_effects <- gasperini_results$effect.list

# create data frame for plotting
plot_df <- data.frame(
  value = c(interaction_effects, enhancer_effects),
  distribution = c(
    rep('Interaction', length(interaction_effects)),
    rep('Enhancer', length(enhancer_effects))
  )
)

# plot boxplot
plot <- ggplot(plot_df, aes(
  x = distribution, 
  y = value
  )) +
  geom_boxplot(fill = '#999999') +
  theme_classic() +
  xlab('') +
  ylab('Effect Size') +
  scale_y_continuous(expand = c(0.02, 0)) +
  theme(
    axis.line = element_line(linewidth = 1),
    axis.title.x = element_text(size = 20, color = 'black'),
    axis.title.y = element_text(size = 20, color = 'black'),
    axis.text = element_text(size = 20, color = 'black'),
    axis.ticks = element_line(color = 'black', linewidth = 1),
    axis.ticks.length = unit(2, 'mm'),
    plot.margin = rep(unit(10, 'mm'), 4),
  )

# save plot to output file
ggsave(
  filename = 'out/boxplot_interaction_enhancer_effects.png',
  device = 'png'
)






