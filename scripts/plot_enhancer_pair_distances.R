# This program plots the distances between enhancers for both the 330 enhancer
# pairs tested and the 3,808 enhancer pairs tested. 
#
# Author: Karthik Guruvayurappan

library(ggplot2)
library(RColorBrewer)

# read in 330 enhancer pair distances
enhancer.distances <- read.csv('/iblm/netapp/data1/external/Gasperini2019/processed/enhancer_distance_330_pairs.csv')

plot <- ggplot(enhancer.distances, aes(x = distance)) +
    geom_histogram(color = 'black') +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    xlab(bquote("Distance (bp)")) + 
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
    filename = '/iblm/netapp/home/karthik/GLiMMIRS/plots/23_04_24_330_enhancer_distance.pdf',
    device = 'pdf',
    plot = plot
)

ggsave(
    filename = '/iblm/netapp/home/karthik/GLiMMIRS/plots/23_04_24_330_enhancer_distance.png',
    device = 'png',
    plot = plot
)

# read in 3,808 enhancer distances
enhancer.distances <- read.csv('/iblm/netapp/data1/external/Gasperini2019/processed/enhancer_distance_at_scale_pairs.csv')

plot <- ggplot(enhancer.distances, aes(x = distance)) +
    geom_histogram(color = 'black') +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    xlab(bquote("Distance (bp)")) + 
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
    filename = '/iblm/netapp/home/karthik/GLiMMIRS/plots/23_04_24_3808_enhancer_distance.pdf',
    device = 'pdf',
    plot = plot
)

ggsave(
    filename = '/iblm/netapp/home/karthik/GLiMMIRS/plots/23_04_24_3808_enhancer_distance.png',
    device = 'png',
    plot = plot
)
