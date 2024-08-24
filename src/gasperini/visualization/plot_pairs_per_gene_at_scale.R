# This script plots the number of enhancer pairs per gene for the at-scale
# analysis.
#
# Author: Karthik Guruvayurappan

library(stats)
library(dplyr)
library(ggplot2)

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

# read in the "true" double perturbation counts
perturbation_counts <- read.csv(
  paste0(
    'data/gasperini/processed/',
    'enhancer_pair_efficiency_adjusted_double_perturb_counts.csv'
  )
)

# merge model outputs with perturbation counts
models <- merge(
  models,
  perturbation_counts,
  by.x = c('enhancer.1.list', 'enhancer.2.list', 'gene.list'),
  by.y = c('enhancer_1_list', 'enhancer_2_list', 'gene_list')
)

# filter for true 10 cell threshold
models <- models[models$double_perturbation_counts >= 10, ]

# count number of pairs per gene
plot_df <- models %>% count(gene.list)

plot <- ggplot(plot_df, aes(x = n)) +
  geom_histogram(fill = 'darkgray', color = 'darkgray') +
  theme_classic() +
  scale_x_continuous(expand = c(0.02, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  xlab("Number of Enhancer Pairs") + 
  ylab("Count") +
  theme(
    axis.line = element_line(linewidth = 1),
    axis.title.x = element_text(size = 24, color = 'black'),
    axis.title.y = element_text(size = 24, color = 'black'),
    axis.text = element_text(size = 24, color = 'black'),
    axis.ticks = element_line(color = 'black', linewidth = 1),
    axis.ticks.length = unit(2, 'mm'),
    plot.margin = rep(unit(10, 'mm'), 4),
  )

ggsave(
  plot = plot,
  filename = 'out/pairs_per_gene_at_scale.pdf',
  device = 'pdf'
)

