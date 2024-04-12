# This script plots histograms of the Cook's distances for cells with the
# double perturbation for the 2 significant interactions observed in the
# STING-seq data.
#
# Author: Karthik Guruvayurappan

library(ggplot2)

# read in first pair Cook's distances
cooks_distances <- read.csv('data/sting_seq/processed/rs1326279_rs1926231_cooks_distances.csv')
cooks_distances <- cooks_distances[cooks_distances$perturbation == 'E1 + E2', ]

plot <- ggplot(cooks_distances, aes(x = cooks_distances)) +
  geom_histogram() +
  theme_classic()

ggsave(
  plot = plot,
  filename = 'out/histogram_sting_seq_rs1326279_rs1926231_cooks_distances.png',
  device = 'png'
)

# read in first pair Cook's distances
cooks_distances <- read.csv('data/sting_seq/processed/rs1926231_rs6669994_cooks_distances.csv')
cooks_distances <- cooks_distances[cooks_distances$perturbation == 'E1 + E2', ]

plot <- ggplot(cooks_distances, aes(x = cooks_distances)) +
  geom_histogram() +
  theme_classic()

ggsave(
  plot = plot,
  filename = 'out/histogram_sting_seq_rs1926231_rs6669994_cooks_distances.png',
  device = 'png'
)
