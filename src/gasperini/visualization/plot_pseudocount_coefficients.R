# Plot the interaction coefficients with and without the use of a pseudocount.
# This uses the 264 pairs analysis.
#
# Author: Karthik Guruvayurappan

library(ggplot2)
library(RColorBrewer)

# read in pseudocount models
pseudocount_models <- read.csv(
  'data/gasperini/processed/enhancer_pairs_330_models.csv'
)

# isolate necessary columns and rename
pseudocount_models <- pseudocount_models[
  ,
  c(
    'enhancer.1.list',
    'enhancer.2.list',
    'gene.list',
    'interaction.effects'
  )
]
colnames(pseudocount_models) <- c(
  'enhancer_1',
  'enhancer_2',
  'gene',
  'interaction_pseudocount'
)

# read in no pseudocount models
no_pseudocount_models <- read.csv(
  'data/gasperini/processed/enhancer_pairs_330_models_no_pseudocount.csv'
)

# isolate necessary columns and rename
no_pseudocount_models <- no_pseudocount_models[
  ,
  c(
    'enhancer.1.list',
    'enhancer.2.list',
    'gene.list',
    'interaction.effects'
  )
]
colnames(no_pseudocount_models) <- c(
  'enhancer_1',
  'enhancer_2',
  'gene',
  'interaction_no_pseudocount'
)

# merge data frames to make plotting df
plot_df <- merge(pseudocount_models, no_pseudocount_models)
plot_df <- plot_df[complete.cases(plot_df), ]

plot <- ggplot(plot_df, aes(
  x = interaction_no_pseudocount, 
  y = interaction_pseudocount
)) +
  geom_abline(slope = 1, intercept = 0, linewidth = 1) +
  geom_point(color = 'darkgray', size = 3) +
  theme_classic() +
  scale_x_continuous(expand = c(0.02, 0)) +
  scale_y_continuous(expand = c(0.02, 0)) +
  xlab('No Pseudocount') +
  ylab('Pseudocount') +
  ggtitle('Interaction Coefficients') +
  theme(
    axis.line = element_line(linewidth = 1),
    axis.title.x = element_text(size = 24, color = 'black'),
    axis.title.y = element_text(size = 24, color = 'black'),
    axis.text = element_text(size = 24, color = 'black', family = 'Helvetica'),
    axis.ticks = element_line(color = 'black', linewidth = 1),
    axis.ticks.length = unit(2, 'mm'),
    plot.margin = rep(unit(10, 'mm'), 4),
    plot.title = element_text(size = 24, color = 'black', hjust = 0.5)
  )

ggsave(
    filename = 'out/scatterplot_pseudocount_coefficients.pdf',
    device = 'pdf',
    plot = plot
)
