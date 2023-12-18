# This script plots the null permutation coefficients for the significant
# STING-seq interactions.
#
# Author: Karthik Guruvayurappan

library(ggplot2)
library(RColorBrewer)
library(gridExtra)

permutation.coefficients <- read.csv('/iblm/netapp/home/karthik/GLiMMIRS/data/experimental/processed/rs1926231_rs6669994_permutation_coefficients.csv')

plot <- ggplot(permutation.coefficients, aes(x = x)) +
    geom_histogram(color = 'black') +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    xlab(bquote("rs1926231 rs6669994")) + 
    ylab(bquote(Count)) +
    theme_classic() +
    theme(
    axis.line = element_line(linewidth = 1),
    axis.title.x = element_text(size = 20, color = 'black'),
    axis.title.y = element_text(size = 20, color = 'black'),
    axis.text = element_text(size = 20, color = 'black'),
    axis.ticks = element_line(color = 'black', linewidth = 1),
    axis.ticks.length = unit(2, 'mm'),
    plot.margin = rep(unit(10, 'mm'), 4),
)

ggsave(
    filename = '/iblm/netapp/home/karthik/GLiMMIRS/out/23_12_11_sting_seq_rs1926231_rs6669994_permutation_coefficients.png',
    device = 'png',
    plot = plot
)

permutation.coefficients <- read.csv('/iblm/netapp/home/karthik/GLiMMIRS/data/experimental/processed/rs1926231_rs1326279_permutation_coefficients.csv')

plot <- ggplot(permutation.coefficients, aes(x = x)) +
    geom_histogram(color = 'black') +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    xlab(bquote("rs1926231 rs1326279")) + 
    ylab(bquote(Count)) +
    theme_classic() +
    theme(
    axis.line = element_line(linewidth = 1),
    axis.title.x = element_text(size = 20, color = 'black'),
    axis.title.y = element_text(size = 20, color = 'black'),
    axis.text = element_text(size = 20, color = 'black'),
    axis.ticks = element_line(color = 'black', linewidth = 1),
    axis.ticks.length = unit(2, 'mm'),
    plot.margin = rep(unit(10, 'mm'), 4),
)

ggsave(
    filename = '/iblm/netapp/home/karthik/GLiMMIRS/out/23_12_11_sting_seq_rs1926231_rs1326279_permutation_coefficients.png',
    device = 'png',
    plot = plot
)
