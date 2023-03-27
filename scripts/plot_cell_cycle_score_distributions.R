# This program plots the distribution of cell cycle scores output by Seurat which will be
# incorporated into the linear model.
#
# This program was written by Karthik Guruvayurappan.

library(ggplot2)
library(RColorBrewer)

s.scores <- read.csv('/iblm/netapp/home/karthik/GLiMMIRS/gasperini_data/s_scores.csv')$S.Score
g2m.scores <- read.csv('/iblm/netapp/home/karthik/GLiMMIRS/gasperini_data/g2m_scores.csv')$G2M.Score
cell.cycle.scores <- data.frame(cbind(s.scores, g2m.scores))

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
    filename = '/iblm/netapp/home/karthik/GLiMMIRS/plots/23_03_26_cell_cycle_s_scores.pdf',
    device = 'pdf',
    plot = plot
)

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
    filename = '/iblm/netapp/home/karthik/GLiMMIRS/plots/23_03_26_cell_cycle_g2m_scores.pdf',
    device = 'pdf',
    plot = plot
)
