# This program plots the distribution of guide efficiencies calculated using the GuideScan 2.0 tool.
# This program was written by Karthik Guruvayurappan.

library(ggplot2)
library(RColorBrewer)

# read in guidescan output
guidescan.guide.info <- read.csv('/iblm/netapp/home/karthik/GuideScan/Gasperini2019/enhancer_guidescan_output.csv')

plot <- ggplot(guidescan.guide.info, aes(x = Cutting.Efficiency)) +
    geom_histogram(color = 'black') +
    scale_x_continuous(expand = c(0, 0), limits = c(0,1)) +
    scale_y_continuous(expand = c(0, 0)) +
    xlab(bquote("Guide Efficiency")) + 
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
    filename = '/iblm/netapp/home/karthik/GLiMMIRS/plots/23_03_25_guide_efficiency_distribution.pdf',
    device = 'pdf',
    plot = plot
)

ggsave(
    filename = '/iblm/netapp/home/karthik/GLiMMIRS/plots/23_03_25_guide_efficiency_distribution.png',
    device = 'png',
    plot = plot
)
