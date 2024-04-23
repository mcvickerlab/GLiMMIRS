# This program plots the distribution of guide efficiencies 
# calculated using the GuideScan 2.0 tool.
#
# This program was written by Karthik Guruvayurappan.

library(stats)
library(ggplot2)
library(RColorBrewer)

# read in guidescan output
guide_info <- read.csv('data/gasperini/interim/enhancer_guide_info.csv')
guide_info <- guide_info[complete.cases(guide_info), ]

plot <- ggplot(guide_info, aes(x = Cutting.Efficiency)) +
    geom_histogram(color = 'black', fill = 'darkgray', binwidth = 0.04) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    xlab(bquote("Guide Efficiency")) + 
    ylab(bquote(Count)) +
    theme_classic() +
    theme(
    axis.line = element_line(linewidth = 0.6),
    axis.title.x = element_text(size = 10, color = 'black'),
    axis.title.y = element_text(size = 10, color = 'black'),
    axis.text = element_text(size = 10, color = 'black'),
    axis.ticks = element_line(color = 'black', linewidth = 0.6),
    axis.ticks.length = unit(2, 'mm'),
    plot.margin = rep(unit(5, 'mm'), 4),
    )
    
ggsave(
    filename = 'out/histogram_guide_efficiency_distribution.pdf',
    device = 'pdf',
    plot = plot,
    unit = 'in',
    width = 3.65,
    height = 2.14
)