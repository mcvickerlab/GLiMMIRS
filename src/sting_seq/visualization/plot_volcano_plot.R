# This script generates a volcano plot showing the interaction effects and
# p-values from running GLiMMIRS on the STING-seq data.
#
# Author: Karthik Guruvayurappan

library(stats)
library(ggplot2)
library(RColorBrewer)

# read in results from STING-seq analysis
results <- read.csv('data/sting_seq/processed/GLiMMIRS_int_ptprc.csv')

# read in Cook's Distances for significant interactions
interaction_1_cooks <- read.csv(
  'data/sting_seq/processed/rs1326279_rs1926231_cooks_distances.csv'
)
interaction_2_cooks <- read.csv(
  'data/sting_seq/processed/rs1926231_rs6669994_cooks_distances.csv'
)

# determine if interaction 1 is filtered
interaction_1_max <- max(interaction_1_cooks$cooks_distances)
interaction_1_mean <- mean(interaction_1_cooks$cooks_distances[
  interaction_1_cooks$perturbation == 'E1 + E2'
])
print(interaction_1_max / interaction_1_mean)
print('Interaction 1 is discarded')

# determine if interaction 2 is filtered
interaction_2_max <- max(interaction_2_cooks$cooks_distances)
interaction_2_mean <- mean(interaction_2_cooks$cooks_distances[
  interaction_2_cooks$perturbation == 'E1 + E2'
])
print(interaction_2_max / interaction_2_mean)
print('Interaction 2 is discarded')

# make labels for plotting
results$labels <- rep('Insignificant', nrow(results))
results$labels[results$fdr_int_pvalues < 0.1] <- 'Outlier'


# plot volcano plot
plot <- ggplot(results, aes(
  x = interaction_estimates, 
  y = -log10(interaction_pvalues),
  color = labels
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
    legend.position = c(0.20, 0.89),
    legend.text = element_text(size = 20)
  ) +
  labs(color = NULL) +
  scale_color_manual(
    breaks = c('Outlier', 'Significant (FDR < 0.1)', 'Insignificant'),
    values = c('deepskyblue', 'red', 'darkgray')
  )

# save plot to output file
ggsave(
  filename = 'out/sting_seq_volcano_plot.pdf',
  plot = plot,
  device = 'pdf',
  width = 7,
  height = 7,
  unit = 'in'
)
