# For a given enhancer-gene pair in the previously published 664 enhancer-gene
# pairs from Gasperini et al., this script plots the difference between the max
# guide efficiency and the minimum guide efficiency and the corresponding
# difference in effect sizes.
#
# Author: Karthik Guruvayurappan

library(dplyr)
library(ggplot2)

# read in results from guide-level models
guide.level.models <- read.csv('/iblm/netapp/data1/external/Gasperini2019/processed/23_12_11_enhancer_gene_pairs_suppl_table_2_baseline_model_guide_level.csv')

# get maximum guide efficiency for each enhancer-gene pair
max.efficiency <- guide.level.models %>%
    group_by(enhancer.list, gene.list) %>%
    filter(efficiency.list == max(efficiency.list, na.rm = TRUE))
max.efficiency <- max.efficiency[, c('enhancer.list', 'gene.list', 'efficiency.list', 'effect.list')]
colnames(max.efficiency) <- c('enhancer', 'gene', 'max.efficiency', 'max.effect')

# get minimum guide efficiency for each enhancer-gene pair
min.efficiency <- guide.level.models %>%
    group_by(enhancer.list, gene.list) %>%
    filter(efficiency.list == min(efficiency.list, na.rm = TRUE))
min.efficiency <- min.efficiency[, c('enhancer.list', 'gene.list', 'efficiency.list', 'effect.list')]
colnames(min.efficiency) <- c('enhancer', 'gene', 'min.efficiency', 'min.effect')

# merge the two dataframes
plot.df <- merge(max.efficiency, min.efficiency)

# get rid of rows where max and min are the same
plot.df <- plot.df[plot.df$max.efficiency != plot.df$min.efficiency, ]

# compute delta guide efficiency and delta effect
plot.df$diff.ge <- plot.df$max.efficiency - plot.df$min.efficiency
plot.df$diff.effect <- plot.df$max.effect - plot.df$min.effect

# plot a scatterplot
plot <- ggplot(plot.df, aes(x = diff.ge, y = diff.effect)) +
    geom_point() +
    theme_classic()

ggsave(
    filename = 'out/23_12_24_delta_guide_efficiency_effect.png',
    device = 'png',
    plot = plot
)

plot.df <- plot.df[abs(plot.df$diff.effect) < 10, ]

plot <- ggplot(plot.df, aes(x = diff.ge, y = diff.effect)) +
    geom_point() +
    theme_classic()

ggsave(
    filename = 'out/23_12_24_delta_guide_efficiency_effect_filtered.png',
    device = 'png',
    plot = plot
)
