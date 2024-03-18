# This script plots a QQ-plot of the interaction term p-values for the 330
# enhancer-enhancer-gene sets where each enhancer comes from the 664
# enhancer-gene pairs previously published by Gasperini et al.
#
# Author: Karthik Guruvayurappan

library(stats)
library(ggplot2)

# read in results from 330 enhancer-enhancer-gene sets
results <- read.csv('data/gasperini/processed/enhancer_pairs_330_models.csv')

# filter out NA values
results <- results[complete.cases(results), ]

# print out dimensions of the results after filtering
print(dim(results))

# sort resulting interaction p-values
results <- results[order(results$interaction.pvalues), ]

# add a uniform distribution of p-values
results$unif <- 1:nrow(results) / nrow(results)

# perform FDR correction on interaction p-values
results$adj.interaction.pvalues <- p.adjust(
  results$interaction.pvalues,
  method = 'fdr'
)

# print number of significant interactions
print(sum(results$adj.interaction.pvalues < 0.1))

# create QQ-plot
plot <- ggplot(results, aes(
  x = -log10(unif), 
  y = -log10(interaction.pvalues)
  )) +
  geom_point(size = 3) +
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
    plot.margin = rep(unit(10, 'mm'), 4)
  )


ggsave(
  filename = 'out/qqplot_interaction_pvalues_330_pairs.png',
  plot = plot,
  device = 'png'
)
