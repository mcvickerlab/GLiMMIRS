# This program plots the distribution of cell cycle scores for the STING-seq
# dataset.
#
# Author: Karthik Guruvayurappan

library(ggplot2)
library(RColorBrewer)

# concatenate cell cycle scores into a dataframe
s.scores <- read.csv('/iblm/netapp/home/karthik/GLiMMIRS/data/experimental/interim/23_11_14_sting_seq_v1_cell_cycle_s_scores.csv')$S.Score
g2m.scores <- read.csv('/iblm/netapp/home/karthik/GLiMMIRS/data/experimental/interim/23_11_14_sting_seq_v1_cell_cycle_g2m_scores.csv')$G2M.Score
cell.cycle.scores <- data.frame(cbind(s.scores, g2m.scores))

# plot S Scores
plot <- ggplot(cell.cycle.scores, aes(x = s.scores)) +
    geom_histogram(color = 'black') +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    xlab(bquote("S Score")) + 
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
    filename = '/iblm/netapp/home/karthik/GLiMMIRS/out/23_11_14_sting_seq_v1_cell_cycle_s_scores.png',
    device = 'png',
    plot = plot
)

# plot G2M Scores
plot <- ggplot(cell.cycle.scores, aes(x = g2m.scores)) +
    geom_histogram(color = 'black') +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    xlab(bquote("G2M Score")) + 
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
    filename = '/iblm/netapp/home/karthik/GLiMMIRS/out/23_11_14_sting_seq_v1_cell_cycle_g2m_scores.png',
    device = 'png',
    plot = plot
)

